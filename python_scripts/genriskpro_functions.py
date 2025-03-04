# -*- coding: utf-8 -*-

import csv
import sys
import pandas as pd
import numpy as np
import os
import gzip
import paramiko
import re
from scp import SCPClient
import time
from itertools import islice
from datetime import datetime
import pickle
print("mudule load test")

def parse_clinvar_vcf(clinvar_vcf_path):
    print(f"Parsing ClinVar VCF from {clinvar_vcf_path}...")
    clinvar_dict = {}      
    with gzip.open(clinvar_vcf_path, 'rt') as file:
        for line in file:
            if line.startswith("#"):
                continue  
            cols = line.strip().split("\t")
            chrom, pos, id_, ref, alt, qual, filter_, info = cols[:8]
            # get needed clinvar info
            clinvar_info = {}
            for key in ["CLNSIG", "CLNDN", "CLNHGVS", "CLNSIGINCL", "CLNVC", "GENEINFO","CLNSIGCONF", "CLNREVSTAT", "CLNDNINCL"]:
                match = re.search(fr"{key}=([^;]+)", info)
                value = None    
                if match:
                    value = match.group(1)
                    ## note this is very important that all csq single value should separate not by |
                    # if isinstance(value, list):
                    #     value = "&".join(value)  # 如果 ClinVar 解析出来是 list，则连接
                    #  **CLNDISDB 里的 `,` 也要替换成 `&`**
                    value = value.replace(",", "&").replace("|", "&")
                clinvar_info[key] = value
            # build the key
            key = (chrom, pos, ref, alt)
            clinvar_dict[key] = clinvar_info
    return clinvar_dict


## For REVEL scores refinement from 0.460&0.460&.&. to 0.460
def extract_decimal_from_string(s):
    if not s or not isinstance(s, str):
        return None
    matches = re.findall(r"\d+\.\d+", s)
    if not matches:
        return None
    # return the first matched value
    return matches[0]


def read_db_file(filepath, encoding="ISO-8859-1", sep="\t", fillna_str="No info", drop_allna_cols=False):
    """help to formally read files"""
    df = pd.read_csv(filepath, sep=sep, encoding=encoding)
    if drop_allna_cols:
        df = df.dropna(axis=1, how='all')
    df = df.replace(np.nan, fillna_str)
    return df

#    use the original "if targetGroup in [xxx] => Basic..., elif => Extended..., else => keep"
def make_final_target_group(tg):
    if tg in [
        'Health predipositions/Disease risk',
        'Carrier-screening',
        'Heriditary-cancer risk syndrome',
        'Newborn-screening'
    ]:
        return 'Basic (for healthy subjects),Usually used for:' + tg
    elif tg in ['Expanded_clinvar_mono_diseases','Expanded_mono_rare_diseases']:
        return 'Extended (for potential patients),' + tg
    else:
        return tg
    
    # generate HGVS standard naming
def generate_hgvs_naming(row):
    hgvsg = row.get('HGVSg', '')  # if there is no HGVSg column, set it to empty
    hgvsc = row.get('HGVSc', '')
    hgvsp = row.get('HGVSp', '')
    gene = row.get('Genes', '')

    # check if HGVSc is not empty
    if pd.notna(hgvsc) and hgvsc.strip():
        try:
            transcript, c_change = hgvsc.split(':', 1)  # extract the transcript and mutation information
            if pd.notna(hgvsp) and hgvsp.strip():
                # rule 1: there is protein change
                p_change = hgvsp.split(':', 1)[1] if ':' in hgvsp else hgvsp
                hgvs_name = f"{transcript}({gene}):{c_change} ({p_change})"
            else:
                # rule 2: there is no protein change
                hgvs_name = f"{transcript}({gene}):{c_change}"
        except ValueError:
            # if the HGVSc format is not correct, skip
            hgvs_name = hgvsg if pd.notna(hgvsg) else ''
    else:
        # rule 3: not in the transcript region
        hgvs_name = hgvsg if pd.notna(hgvsg) else ''

    return hgvs_name


# generate ReviewStar column
def add_review_star(df):
    # generate ReviewStar based on ClinVar_CLNREVSTAT
    def calculate_review_star(review):
    # if it is empty (None/NaN), convert it to string
        review = "" if (review is None or pd.isna(review)) else str(review)
        r = review.lower().replace(" ", "").replace("_", "").replace("&", "").replace(",", "")
        if r in ("criteriaprovidedconflictinginterpretations",
                "criteriaprovidedconflictingclassifications"):
            return 1
        elif r == "criteriaprovidedsinglesubmitter":
            return 1
        elif r == "criteriaprovidedmultiplesubmittersnoconflicts":
            return 2
        elif r == "reviewedbyexpertpanel":
            return 3
        elif r == "practiceguideline":
            return 4
        else:
            return 0
    # overwrite the existing ReviewStar column
    df['ReviewStar'] = df['ClinVar_CLNREVSTAT'].apply(calculate_review_star)
    return df

# 过滤Pathogenic和Likely_pathogenic变异

def filter_known_benign_variants(df):
    # ClinVar_CLNSIGCONF needs to include "Pathogenic" or "Likely_pathogenic"，but can not include "benign"
    pathogenic_mask = df['ClinVar_CLNSIGCONF'].str.contains('Pathogenic|Likely_pathogenic', na=False, case=False) & \
                      ~df['ClinVar_CLNSIGCONF'].str.contains('benign', na=False, case=False)

    # if ClinVar_CLNSIGCONF is empty, check CLIN_SIG if include "benign" or "likely_benign" and remove
    clin_sig_mask = ~df['CLIN_SIG'].str.contains('benign|likely_benign', na=False, case=False)

    # merge two criteria 
    # 1. if ClinVar_CLNSIGCONF ok，keep
    # 2. if ClinVar_CLNSIGCONF empty，check CLIN_SIG，to exclude benign
    filtered_df = df[pathogenic_mask | (df['ClinVar_CLNSIGCONF'].isna() & clin_sig_mask)]
    
    return filtered_df

def get_apoe_isoform_and_risk(rs429358_genotype, rs7412_genotype):
    """
    based on rs429358_genotype and rs7412_genotype
    return (APOE Isoform, AD risk description) tuple.

    reference table:
    rs429358(Ref:T)  rs7412(Ref:C)   Isoform   AD risk
**  TT        TT       ε2/ε2    Reduced Risk (protective)
?   TT        CT       ε2/ε3    Neutral to Reduced Risk
?   CC        CT       ε2/ε3    Neutral to Reduced Risk
**  TT        CC       ε3/ε3    Neutral (most common)  
?   TC        CC       ε3/ε4    Increased Risk
**  CC        CC       ε4/ε4    Highest Risk
    """
    genotype_map = {
        ("TT", "TT"): ("e2/e2", "Reduced Risk (protective)"),
        ("TT", "CT"): ("e2/e3", "Neutral to Reduced Risk"),
        ("TT", "CC"): ("e3/e3", "Neutral (most common)"),
        ("TC", "CC"): ("e3/e4", "Increased Risk"),
        ("CC", "CC"): ("e4/e4", "Highest Risk"),
    }
    return genotype_map.get((rs429358_genotype, rs7412_genotype), ("Unknown", "Unknown"))

def parse_genotype(entry: str, ref: str, alt: str, chrom: str, check_male: bool) -> tuple:
    """
    Parse genotype information and return genotype and zygosity
    
    Args:
        entry: Sample field string (e.g. "0/1:25,6:31...")
        ref: REF field value
        alt: ALT field value 
        chrom: Chromosome field value
        check_male: Flag to check male samples
        
    Returns:
        tuple: (genotype string, zygosity)
    """
    # initialize default values
    genotype = "Missing"
    zygosity = "Unknown"
    
    # handle empty value case
    if not entry or entry == './.':
        return (genotype, zygosity)
    
    try:
        # split the genotype part
        gt_part = entry.split(':', 1)[0]
        if '.' in gt_part:
            return ("Missing", "Unknown")
        
        # parse the allele list
        ref_alleles = ref.split(',')
        alt_alleles = alt.split(',')
        all_alleles = ref_alleles + alt_alleles
        
        # parse the genotype index
        idx1, idx2 = map(int, gt_part.replace('/', '|').split('|')[:2])
        
        # build the genotype description
        base_allele = ref_alleles[0]
        observed = []
        for idx in [idx1, idx2]:
            if 0 <= idx < len(all_alleles):
                observed.append(all_alleles[idx])
            else:
                observed.append('N')    # represent invalid index
        
        # generate the genotype string
        original = f"{base_allele}/{base_allele}"
        observed_str = '/'.join(observed)
        genotype = f"{original} > {observed_str}"
        
        # determine the zygosity
        if idx1 == idx2:
            zygosity = 'Homozygous'
        elif idx1 == 0 or idx2 == 0:
            zygosity = 'Heterozygous'
        elif idx1 >= len(ref_alleles) or idx2 >= len(ref_alleles):
            zygosity = 'Compound heterozygous'
        
        # handle the special case of sex chromosome
        if check_male and chrom in ('chrX', 'X') and zygosity == 'Heterozygous':
            zygosity = 'Hemizygous'
            
    except (ValueError, IndexError, AttributeError) as e:
        print(f"Error parsing genotype: {str(e)}")
        return ("Invalid", "Error")
    
    return (genotype, zygosity)