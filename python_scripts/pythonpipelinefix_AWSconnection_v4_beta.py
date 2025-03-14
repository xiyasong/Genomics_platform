#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:41:09 2023

@author: xiyas
"""

#%%Cell 1 Define the path for parameters and databases
# ==============================================================================
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
import genriskpro_functions as gpf

all_start_time = time.time()
# ====== 1) setting up paths and running parameters ======
path = 'local'  # local / uppmax / SZA / sysmed

run_GWAS = False
run_Pharmaco = False
run_traits = False
Upload_To_Cloud = False
Update_Clinvar = False

# config file of setting paths
config = {
    "local": {
        "fileName": "/Users/xiyas/Downloads/bc_40_madeleine_merged2.annotated.vcf.gz",
        "outFile": "/Users/xiyas/ATPM_Project/ML_scores/clinvar_20240611.annotated.vcf.gz.txt",
        "geneBaseFile": "/Users/xiyas/ATPM_Project/database-v3/GeneDB_GenCC.txt",
        "diseaesDBfile": "/Users/xiyas/V2_Genome_reporting/database-file-v2/diseaseDB_1115_3.txt",
        "OMIM_inheritance_DBfile": "/Users/xiyas/V2_Genome_reporting/database-file-v2/pheno_OMIM_all.txt",
        "clingen_file": "/Users/xiyas/ATPM_Project/database-v3/Clingen-variants-2024-12-09.txt",
        "ontology_file": "/Users/xiyas/ATPM_Project/database-v3/genedb.ontology.all0307.csv",
        "GWAS_dbfile": "/Users/xiyas/ATPM_Project/database-v3/Merged_GWAS_vcf_2024.txt",
        "Pharma_dbfile": "/Users/xiyas/ATPM_Project/database-v3/Merged_Pharma_vcf_2024.txt",
        "Trait_dbfile": "/Users/xiyas/ATPM_Project/database-v3/Reports_genome_databases_traits_merged_2.txt",
        "clinvar_vcf_path": "/Users/xiyas/ATPM_Project/database-v3/clinvar_20240611.vcf.gz"
    },
    "uppmax": {
        "fileName": None,  # will read from sys.argv[1] 
        "outFile": None,   # will read from sys.argv[2] 
        "geneBaseFile": "/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/GeneDB.txt",
        "diseaesDBfile": "/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/diseaseDB_1115_3.txt",
        "OMIM_inheritance_DBfile": "/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/pheno_OMIM_all.txt",
        "GWAS_dbfile": "/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Merged_GWAS_vcf.txt",
        "Pharma_dbfile": "/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Merged_Pharma_vcf.txt",
        "Trait_dbfile": "/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Reports_genome_databases_traits_merged.txt",
    },
    "SZA": {
        "fileName": None,
        "outFile": None,
        "geneBaseFile": "/mnt/SZAPORTAL/database-file-v2/GeneDB.txt",
        "diseaesDBfile": "/mnt/SZAPORTAL/database-file-v2/diseaseDB_1115_3.txt",
        "OMIM_inheritance_DBfile": "/mnt/SZAPORTAL/database-file-v2/pheno_OMIM_all.txt",
        "GWAS_dbfile": "/sza_data/GenomicPipeline/database-v3/Merged_GWAS_vcf_2024.txt",
        "Pharma_dbfile": "/sza_data/GenomicPipeline/database-v3/Merged_Pharma_vcf_2024.txt",
        "Trait_dbfile": "/sza_data/GenomicPipeline/database-v3/Reports_genome_databases_traits_merged.txt",
    },
    "sysmed": {
        "fileName": None,
        "outFile": None,
        "clingen_file": "/mnt/storage_pool/Genomics/Genome/database-files/Clingen-variants-2024-12-09.txt",
        "geneBaseFile": "/mnt/storage_pool/Genomics/Genome/database-files/GeneDB_GenCC.txt",
        "diseaesDBfile": "/mnt/storage_pool/Genomics/Genome/database-files/diseaseDB_1115_3.txt",
        "OMIM_inheritance_DBfile": "/mnt/storage_pool/Genomics/Genome/database-files/pheno_OMIM_all.txt",
        "ontology_file": "/mnt/storage_pool/Genomics/Genome/database-files/genedb.ontology.all0307.csv",
        "GWAS_dbfile": "/mnt/storage_pool/Genomics/Genome/database-files/Merged_GWAS_vcf_2024.txt",
        "Pharma_dbfile": "/mnt/storage_pool/Genomics/Genome/database-files/Merged_Pharma_vcf_2024.txt",
        "Trait_dbfile": "/mnt/storage_pool/Genomics/Genome/database-files/Reports_genome_databases_traits_merged_2.txt",
        "clinvar_vcf_path": "/mnt/storage_pool/Genomics/Genome/database-files/clinvar_20240611.vcf.gz"
    }
}

# ====== 2) Select path and loading ======
if path not in config:
    raise ValueError(f"unknown path: {path}")

cfg = config[path]

#   if fileName / outFile  None，read from sys.argv
if cfg["fileName"] is None:
    # 确保 sys.argv 长度足够
    if len(sys.argv) < 3:
        raise ValueError(f"provide inputFile and outputFile,eg. python script.py input.vcf output.txt")
    fileName = sys.argv[1]
    outFile = sys.argv[2]
else:
    fileName = cfg["fileName"]
    outFile = cfg["outFile"]  
    if outFile is None:
        outFile = "test_output.txt"

print(f"Set path = {path}")
print("fileName:", fileName)
print("outFile:", outFile)

# ====== 3) Reading necessary files ======
DiseaseDB   = gpf.read_db_file(cfg["diseaesDBfile"])
OMIM_Inheritance_DB      = gpf.read_db_file(cfg["OMIM_inheritance_DBfile"])
OMIM_Inheritance_DB['phenotypeMimNumber'] = OMIM_Inheritance_DB['phenotypeMimNumber'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")
geneBaseFile = cfg["geneBaseFile"]
diseaesDBfile = cfg["diseaesDBfile"]

# if clingen_file / ontology_file does not exist in some environment, skip or write as cfg.get("clingen_file", None)
clingen_file = cfg.get("clingen_file")
if clingen_file is not None:
    clingene = pd.read_csv(clingen_file, sep="\t",header=0)
    pass

ontology_file = cfg.get("ontology_file")
if ontology_file is not None:
    ontology = pd.read_csv(ontology_file, sep=",",header = 0)
    pass

# ====== 4) Decision of whether running GWAS/Pharmaco/Trait ======
if run_GWAS:
    GWAS_db = gpf.read_db_file(cfg["GWAS_dbfile"])
    print("GWAS matching will run")

if run_Pharmaco:
    Pharma_db = gpf.read_db_file(cfg["Pharma_dbfile"])
    print("Pharmaco matching will run")

if run_traits:
    Trait_db = gpf.read_db_file(cfg["Trait_dbfile"], drop_allna_cols=True)
    trait_list = Trait_db['variants'].to_list()
    print("Traits matching will run")


# ClinVar DB dynamic updates: parse into clinvar_cache

if Update_Clinvar == True:
    # read and create dict
    clinvar_vcf = cfg.get("clinvar_vcf_path")
    clinvar_cache = cfg.get("clinvar_cache_path")
    clinvar_dict = gpf.parse_clinvar_vcf(clinvar_vcf)


#%%Cell 2 start the analysis program for the vep annotated vcf files, generated report A B C and D || E for the traits
# ==============================================================================
start_time = time.time()

#####check it's male or female, to get the correct Zygosity for the chrX variants
#check_male_flag = 0>1
check_male_flag = False
PATHOGENIC_TERMS = {'pathogenic', 'likely_pathogenic'}
BENIGN_TERMS = {'benign', 'likely_benign'}
#tLine = tLine.decode("ISO-8859-1")
i = 0
reportA ,reportB,reportC,reportD,reportE = [], [], [], [], []

# col_map is used to build a mapping once colNames_CSQ is obtained
col_map = {}  # Will be populated after parsing the "ID=CSQ" line
headings = []

# Reading vcf.gz file 
file = gzip.open(fileName,'rt')
tLine = file.readline()

while tLine:
    # remove the newline character
    tLine = tLine.rstrip('\n')
    # split the current line
    iContent = tLine.split('\t')
    i += 1
    ##get the content from VCF annotation header
    if tLine.startswith('#'):
        if 'ID=CSQ' in tLine:
            annoText = iContent[0].split('Format: ')
            colNames_CSQ = annoText[1].replace('">','')
            colNames_CSQ = colNames_CSQ.split('|')
            # construct col_map for all use
            col_map = { name: idx for idx, name in enumerate(colNames_CSQ) }
        elif tLine.startswith('#CHROM'):
            headings = iContent
        # directly goes into next line
        tLine = file.readline()
        continue
    #start processing real data rows
    if not headings:
        tLine = file.readline()
        continue
    # iContent is a VCF line => [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sampleGenotype...]
    # Get multiple transcripts annotation (iText)
    # CSQ includes in INFO field
    # this is going to analyse all transcript result inside CSQ and any one satisfied the criteria, the whole variant goes into report A/B.
    #### =================ClinVar version updates ##### =================
    if Update_Clinvar == True:
        chrom, pos, id_, ref, alt, qual, filter_, info = iContent[:8]
        info_parts = info.split(";")
        csq_index = next((i for i, part in enumerate(info_parts) if part.startswith("CSQ=")), None)
        csq_data = info_parts[csq_index].replace("CSQ=", "").split(",")
        key = (chrom, pos, ref, alt)
        # **if ClinVar has updates, replace the ClinVar info in `CSQ`**
        if key in clinvar_dict:
            new_clinvar_info = clinvar_dict[key]  # ClinVar latest data
            updated_csq_data = []
            for csq_entry in csq_data:
                csq_fields = csq_entry.split("|")  # split CSQ info
                # loop through ClinVar related fields, and use `col_map` to find the correct index to replace
                for field in clinvar_dict[key]:
                    csq_field_name = field  # use the original field name by default
                    # if the field in `CSQ` is start with `ClinVar_`, automatically adjust
                    if f"ClinVar_{field}" in col_map:
                        csq_field_name = f"ClinVar_{field}"
                    if csq_field_name in col_map and new_clinvar_info.get(field):
                        csq_fields[col_map[csq_field_name]] = new_clinvar_info[field]

                updated_csq_data.append("|".join(csq_fields))

            # rejoin `CSQ` and update `INFO`
            info_parts[csq_index] = "CSQ=" + ",".join(updated_csq_data)
            iContent[headings.index('INFO')] = ";".join(info_parts)  # update iContent[7] (INFO)
        #**rejoin `iContent` to `tLine`**
        tLine = "\t".join(iContent)
    #### =================ClinVar version updates ##### =================
    saveFlag1, saveFlag2,saveFlag3,saveFlag4, saveFlag5 = False, False, False, False, False
    iText = [s for s in iContent[headings.index('INFO')].split(';') if 'CSQ=' in s]
    iText = iText[0].replace('CSQ=','').split(',')
    for j in range(0,len(iText)):
        jText = iText[j].split('|')
        # fixing REVEL score 
        if 'REVEL' in col_map and 'REVEL_score' in col_map:
            revel_idx = col_map['REVEL']
            revel_score_idx = col_map['REVEL_score']
            if jText[revel_idx] == '':
                parsed_value = gpf.extract_decimal_from_string(jText[revel_score_idx])
                if parsed_value is not None:
                    jText[revel_idx] = parsed_value

        # 1) ClinP_LP_var, GWAS_var and Pharma_var
        if 'Database_SZAID' in col_map and 'Database_Type' in col_map:
            db_szaid_val = jText[col_map['Database_SZAID']]
            db_type_val  = jText[col_map['Database_Type']]
            if 'ClinVar_CLNSIG' in col_map:
                jCLNSIG = jText[col_map['ClinVar_CLNSIG']]
            # save known pathogenicrecords to reportA, and also it is based on the clinvar_vcf——20240611 (updated)
            if 'Conflicting' in jCLNSIG:
                raw_pathogenicity = jText[col_map['ClinVar_CLNSIGCONF']].split('&')
                pathogenicity = [
                    re.sub(r'\(\d+\)', '', s).strip().lower()
                    for s in raw_pathogenicity
                ]
                has_pathogenic = any(p in PATHOGENIC_TERMS for p in pathogenicity)
                no_benign = not any(b in BENIGN_TERMS for b in pathogenicity)
                CPLP_condition = has_pathogenic and no_benign
            else:
                CPLP_condition = False
            if ('ClinP_LP_var' in db_type_val.split('&') or 
                any(p in jCLNSIG.lower() for p in PATHOGENIC_TERMS) or 
                CPLP_condition):
                saveFlag1 = 'ClinP_LP_var'

        if run_GWAS and 'Database_Type' in col_map:
            if 'GWAS_var' in db_type_val.split('&'):
                saveFlag3 = 'GWAS_var'

        if run_Pharmaco and 'Database_Type' in col_map:
            if 'Pharma_var' in db_type_val.split('&'):
                saveFlag4 = 'Pharma_var'

        # 2) putative/predicted risk variants : not saveFlag1 + MAX_AF<0.05 + 其他条件 （update here:remove Benign/Likely_benign）
        if not saveFlag1 and 'MAX_AF' in col_map:
            max_af_val = jText[col_map['MAX_AF']]
            # clinvar filter
            clinvar_sig = ''
            if 'ClinVar_CLNSIG' in col_map:
                clinvar_sig = jText[col_map['ClinVar_CLNSIG']].lower()  
            # 复合条件判断（增加排除Benign/Likely_benign）
            if (max_af_val == '' or float(max_af_val) < 0.05) and \
               not re.search(r'benign|likely_benign', clinvar_sig):
                # IMPACT / ada_score / rf_score / REVEL / SpliceAI_pred... / BayesDel... / am_class / am_pathogenicity / LoF
                predicted_impact = (
                    ('IMPACT' in col_map and jText[col_map['IMPACT']] == 'HIGH')
                    or ('ada_score' in col_map and jText[col_map['ada_score']] != '' and float(jText[col_map['ada_score']]) > 0.6)
                    or ('rf_score' in col_map and jText[col_map['rf_score']] != '' and float(jText[col_map['rf_score']]) > 0.6)
                    or ('REVEL' in col_map and jText[col_map['REVEL']] != '' and float(jText[col_map['REVEL']]) > 0.75)
                    or ('SpliceAI_pred_DS_AL' in col_map and jText[col_map['SpliceAI_pred_DS_AL']] != '' and float(jText[col_map['SpliceAI_pred_DS_AL']])>0.5)
                    or ('SpliceAI_pred_DS_DG' in col_map and jText[col_map['SpliceAI_pred_DS_DG']] != '' and float(jText[col_map['SpliceAI_pred_DS_DG']])>0.5)
                    or ('SpliceAI_pred_DS_DL' in col_map and jText[col_map['SpliceAI_pred_DS_DL']] != '' and float(jText[col_map['SpliceAI_pred_DS_DL']])>0.5)
                    or ('SpliceAI_pred_DS_AG' in col_map and jText[col_map['SpliceAI_pred_DS_AG']] != '' and float(jText[col_map['SpliceAI_pred_DS_AG']])>0.5)
                    or ('BayesDel_addAF_score' in col_map and jText[col_map['BayesDel_addAF_score']] != '' and float(jText[col_map['BayesDel_addAF_score']])>0.0692655)
                    or ('BayesDel_noAF_score' in col_map and jText[col_map['BayesDel_noAF_score']] != '' and float(jText[col_map['BayesDel_noAF_score']])>-0.0570105)
                    or ('am_class' in col_map and jText[col_map['am_class']] == 'likely_pathogenic'
                        and 'am_pathogenicity' in col_map
                        and jText[col_map['am_pathogenicity']] != ''
                        and float(jText[col_map['am_pathogenicity']])>0.564)
                    or ('LoF' in col_map and jText[col_map['LoF']] == 'HC')
                )
                if predicted_impact:
                    saveFlag2 = "Putative_var"

        # 3) judgement for Trait, need run_traits=True and trait_list
        if run_traits and 'Existing_variation' in col_map:
            ex_variation_val = jText[col_map['Existing_variation']]
            # maybe "rsXXXXX&COSVxxxxx", first split
            ex_variants = ex_variation_val.split('&')
            if ex_variants and ex_variants[0] in trait_list:
                saveFlag5 = "Trait_var"
    # after for j in range(len(iText)) loop, if saveFlag1/2/3/4/5 has value, append the line to respective report
    if saveFlag1:
        reportA.append(tLine)
    if saveFlag2:
        reportB.append(tLine)
    if saveFlag3:
        reportC.append(tLine)
    if saveFlag4:
        reportD.append(tLine)
    if saveFlag5:
        reportE.append(tLine)

    # detect chrY => gender
    if "chrY" in tLine:
        check_male_flag = "Male"

    # print progress every 1000000 lines
    if i % 1000000 == 0:
        print(f"{i} lines processed!")

    # read the next line ---------------
    tLine = file.readline()

gender = "Male" if check_male_flag else "Female"
if check_male_flag :
    print("Male")
else:
    print("Female")

file.close()
print('Cell 2 VEP annotated File processing done! Now start to map GeneDB and DiseaseDB')
end_time = time.time()
print("Total processing time: {:.2f} seconds".format(end_time - start_time))

#%%Cell 3 read gene db and disease db
# ==============================================================================
print("Cell 3: Reading GeneDB and DiseaseDB")
# read gene based database
ConfidenceLevel = []
TargetGroup = []
geneBasedRef = []
Disease = []
SZAdiseaseID_GDB = []
#inh = []
#mark = []
#mim = []
with open(geneBaseFile,'r') as f:
    for line in f:
        line = line.replace('\n','')
        temp = line.split('\t')
        ConfidenceLevel.append(temp[0])
        TargetGroup.append(temp[1])
        geneBasedRef.append(temp[2])
        Disease.append(temp[3])
        SZAdiseaseID_GDB.append(temp[4])
#        inh.append(temp[5])
#        mark.append(temp[6])
#        mim.append(temp[7])
geneBasedRefHeadings = [ConfidenceLevel[0],TargetGroup[0],geneBasedRef[0],Disease[0]]
geneBasedRef.pop(0)
Disease.pop(0)
ConfidenceLevel.pop(0)
TargetGroup.pop(0)
SZAdiseaseID_GDB.pop(0)
#inh.pop(0)
#mim.pop(0)
#mark.pop(0)

# read disease based database
DiseaseID_DSDB = []
DiseaseName_DSDB = []
count=0
with open(diseaesDBfile,'r',encoding = "ISO-8859-1") as f:
    for line in f:
        count+=1
        line = line.replace('\n','')
        temp = line.split('\t')
        DiseaseID_DSDB.append(temp[0])
        DiseaseName_DSDB.append(temp[1])
        #print(line)
DiseaseID_DSDB.pop(0)
DiseaseName_DSDB.pop(0)
disease_lookup = {}
with open(diseaesDBfile,'r',encoding = "ISO-8859-1") as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # 跳过标题行
    for idx, row in enumerate(reader):
        if len(row) < 2:
            continue
        did = row[0]
        name = row[1]
        # keep the original list
        DiseaseID_DSDB.append(did)
        DiseaseName_DSDB.append(name)
        # build the lookup dictionary (single key)
        standardized = name.replace(',', '').strip().lower()
        disease_lookup[standardized] = (did, name, idx)
        
newContentHeadings = ['SZAID','Database','Database_type','SZAdiseaseID','ClinVar_CLNSIG','ClinVar_CLNDN','SZAreportCategory','ClinVar_ID','selCSQ','Zygosity','Genotype']+geneBasedRefHeadings+headings;
newContent = '\t'.join(newContentHeadings)+'\n'

#%%Cell 4 report A : ClinVar Variants (Known Pathogenicity) =================================================
print("Cell 4 report A is for Clinvar Variants")
content_buffer = []

COLUMN_INDEX = {
    'INFO': headings.index('INFO'),
    'CHROM': headings.index('#CHROM'),
    'REF': headings.index('REF'),
    'ALT': headings.index('ALT'),
    'SAMPLE': len(headings) - 1 if len(headings) > 9 else None
}

for i in range(0,len(reportA)):
    if i % 100 == 0:
        print(f"{i} lines processed!")
    iTemp = reportA[i].replace('\n','').split('\t')
    iTemp[headings.index('INFO')] = iTemp[COLUMN_INDEX['INFO']].split(';CSQ=')[0]
    iText = [s for s in reportA[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
    #iText sepatated different transcript annotations (current vep version only pick one transcripts)

    ### ============================== Genotype processing part: if sample, get iGenotype and iZygo ===============
    iREF = reportA[i].split('\t')[COLUMN_INDEX['REF']].split(',')
    iALT = reportA[i].split('\t')[COLUMN_INDEX['ALT']].split(',')
    iGenoTypeList = iREF+iALT
    
    # 1) check there is no missingness of sample field, if no sample field, skip the genotype processing
    # 2) check there is no missingness in the genotype field, if genotype field is missing, skip the genotype processing
    if COLUMN_INDEX['SAMPLE'] is not None:

        sample_entry = reportA[i].split('\t')[COLUMN_INDEX['SAMPLE']]
        chrom = reportA[i].split('\t')[COLUMN_INDEX['CHROM']]
        ref = reportA[i].split('\t')[COLUMN_INDEX['REF']]
        alt = reportA[i].split('\t')[COLUMN_INDEX['ALT']]

        iGenotype, iZygo = gpf.parse_genotype(
            entry=sample_entry,
            ref=ref,
            alt=alt,
            chrom=chrom,
            check_male=check_male_flag
        )
    else:
        iGenotype = "Missing"
        iZygo = "Unknown"  
    ### ==============================
    for j in range(0,len(iText)):
        jText = iText[j].split('|')
        jGenes1 = jText[3]
        jGenes2Temp = jText[col_map["ClinVar_GENEINFO"]].split('&')
        jGenes2 = []
        for k in range(0,len(jGenes2Temp)):
            ### split by ; or : for ClinVar_GENEINFO new format
            kTemp = re.split('[;:]', jGenes2Temp[k])
            jGenes2.append(kTemp[0])
        jGenes = list(set([jGenes1]+jGenes2))
        jGenes = [gene for gene in jGenes if gene != ""]
        ######### clear GWAS and pharmaco id for clinvar variants ######
        ###The reason to do this: GWAS,pharmaco, Clinvar can be overlapped, so several SZAvar can be actually one, just showing in both clinvar and gwas or ..
        ###This will have issues like ValueError: 'ClinP_LP_var' is not in list
        ###Is because in one of the transcripts, the way of express alternative allele is different with what clinvar used.
        ## expression of "-" will annotated Clinvar_P_LP, while "TGCTGC" will be annotated Clinvar_LP_var
        if jText[col_map['Database_Type']]!= '':
            #if jText[colNames_CSQ.index('Database_Type')] == "GWAS_var&ClinP_LP_var":
            #if jText[colNames_CSQ.index('Database_Type')] == "ClinP_LP_var&GWAS_var":
            jCLINVARind = jText[col_map['Database_Type']].split('&').index('ClinP_LP_var')
            jDatabaseType = jText[col_map['Database_Type']].split('&')[jCLINVARind]
            #jDatabase = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
            jCLNID = jText[col_map['Database']].split('&')[jCLINVARind]
            jSZAvarID = jText[col_map['Database_SZAID']].split('&')[jCLINVARind]
            jCLNSIG = jText[col_map['ClinVar_CLNSIG']]
            #jCLNSIG = jText[colNames_CSQ.index('ClinVar_CLNSIG')].split('&')
            jClinVar_CLNDN =  jText[col_map['ClinVar_CLNDN']].replace('&_',',_').replace('&','|').replace('_',' ')
            if 'Pathogenic' in jCLNSIG:
                scoreFlag = 5
            elif 'Likely_pathogenic' in jCLNSIG:
                scoreFlag = 4
            elif 'Conflicting_classifications_of_pathogenicity' in jCLNSIG or 'Conflicting_interpretations_of_pathogenicity' in jCLNSIG:
                pathogenicity = jText[col_map['ClinVar_CLNSIGCONF']].split('&')
                pathogenicity = [re.sub(r'\(\d+\)', '', s).lower() for s in pathogenicity]
                if ('pathogenic' in pathogenicity) and not any(s in pathogenicity for s in ['benign', 'likely_benign']):
                    scoreFlag = 2
                elif ('likely_pathogenic' in pathogenicity) and not any(s in pathogenicity for s in ['benign', 'likely_benign']):
                    scoreFlag = 1
                else:
                    scoreFlag = 0   
            else:
                scoreFlag = 0
            # Case one: if the associated gene can be mapped to GeneDB
            if scoreFlag > 0 and any([x for x in jGenes if x in geneBasedRef]):
                jGeneList = [x for x in jGenes if x in geneBasedRef]
                jInds = []
                for k in range(0,len(jGeneList)):
                    kTemp = jGeneList[k]
                    jInds = jInds +[ii for ii, x in enumerate(geneBasedRef) if x == kTemp]
                ##### =================================1.Match Disease Name between Clinvar and DiseaseDB, for adding disease info
                kClinVarDiseases = jClinVar_CLNDN.split('|')
                #print(kClinVarDiseases)
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                for k in range(0,len(kClinVarDiseases)):
                    kClinVarDisease_temp = kClinVarDiseases[k].replace(",", "").strip().lower()
                    matched = disease_lookup.get(kClinVarDisease_temp, None)
                    if matched:
                        did, dname, _ = matched
                        jVarSZAdiseaseIDs.append(did)
                        jVarDiseaseNames.append(dname)
                    else:
                        jVarSZAdiseaseIDs.append('')
                        jVarDiseaseNames.append(kClinVarDiseases[k])                   
                    # =================================
                ##### 2. Match Disease Name between Clinvar and GeneDB, for getting reporting confidence score
                ## To give the score on gene level for ClinVar Variants (Known Pathogenicity), should match Clinvar disease type
                ## fixed by if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs: SZAdiseaseID_GDB is the disease from GeneDB, to see
                # whether it's matching with the ClinVar_CLNDN disease jVarSZAdiseaseIDs
                for k in range(0,len(jInds)):
                    if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs:
                        #SZAscore =0
                        if ConfidenceLevel[jInds[k]] == 'High_confidence':
                            #print('SZAdiseaseID_GDB[jInds[k]]=',SZAdiseaseID_GDB[jInds[k]])
                            #print(' jVarSZAdiseaseIDs=',jVarSZAdiseaseIDs)
                            SZAscore = 10+scoreFlag
                        elif ConfidenceLevel[jInds[k]] == 'Moderate_confidence':
                            SZAscore = 5+scoreFlag
                        elif ConfidenceLevel[jInds[k]] == 'Low_confidence':
                            SZAscore = scoreFlag
                        else:
                            sys.exit('Unknown confidence score level for the report!')
                        jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                        row = '\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],Disease[jInds[k]]]+iTemp)+'\n'
                        content_buffer.append(row)
                    else:
                        GeneDB_szaID = SZAdiseaseID_GDB[jInds[k]]
                        GeneDB_Disease = Disease[jInds[k]]
                        GeneDB_Disease_word = set(GeneDB_Disease.upper().split())

                        for name in range(len(jVarDiseaseNames)):
                            clinvar_disease_word = set(jVarDiseaseNames[name].upper().split())
                            matches = False
                            if clinvar_disease_word.issubset(GeneDB_Disease_word):
                                matches = True
                            elif clinvar_disease_word.issubset(GeneDB_Disease_word):
                                matches = True
                            if matches:
                                if ConfidenceLevel[jInds[k]] == 'High_confidence':
                                #print('SZAdiseaseID_GDB[jInds[k]]=',SZAdiseaseID_GDB[jInds[k]])
                                #print(' jVarSZAdiseaseIDs=',jVarSZAdiseaseIDs)
                                    SZAscore = 10 + scoreFlag
                                elif ConfidenceLevel[jInds[k]] == 'Moderate_confidence':
                                    SZAscore = 5 + scoreFlag
                                elif ConfidenceLevel[jInds[k]] == 'Low_confidence':
                                    SZAscore = scoreFlag
                                else:
                                    sys.exit('Unknown confidence score level for the report!')
                                
                                jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                                row = '\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                                Disease[jInds[k]]]+iTemp)+'\n'
                                content_buffer.append(row)

                        else:
                            ########here, if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs is false, then no SZAscore value
                            SZAscore = scoreFlag
                            jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                            row = '\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                            Disease[jInds[k]]]+iTemp)+'\n'
                            content_buffer.append(row)
            # Case two: if there is jGenes, but not included in GeneDB
            elif scoreFlag >0 and any(jGenes): 
                #print("here")
                #sys.exit('Found a gene set with '+', '.join(jGenes)+' that is not included in GeneDB ')
                SZAscore = scoreFlag
                kClinVarDiseases = jClinVar_CLNDN.split('|')
                #print(kClinVarDiseases)
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                jSZADiseaseID = []
                jDisease = []
                for k in range(0,len(kClinVarDiseases)):
                    kClinVarDisease_temp = kClinVarDiseases[k].replace(",", "").strip().lower()
                    matched = disease_lookup.get(kClinVarDisease_temp, None)
                    if matched:
                        did, dname, _ = matched
                        jVarSZAdiseaseIDs.append(did)
                        jVarDiseaseNames.append(dname)
                    else:
                        jSZADiseaseID.append('')  # 追加空字符串作为占位符
                        jDisease.append('')  
                ### handle multiple diseases matching
                jDisease = [x for x in jDisease if x != '']
                jSZADiseaseID = [x for x in jSZADiseaseID if x != '']
                if len(jDisease) == 1:
                    jDisease = ''.join(jDisease)
                    jSZADiseaseID = ''.join(jSZADiseaseID)
                    row = '\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease]+iTemp)+'\n'
                    content_buffer.append(row)
                elif len(jDisease) > 1:
                    for m in range(len(jDisease)):
                        mjDisease = ''.join(jDisease[m])
                        mjSZADiseaseID = ''.join(jSZADiseaseID[m])
                        row = '\t'.join([jSZAvarID,jCLNID,jDatabaseType,mjSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),mjDisease]+iTemp)+'\n'
                        content_buffer.append(row)         
            # Case three: if there is no jGenes
            elif scoreFlag >0:
                SZAscore = scoreFlag
                kClinVarDiseases = jClinVar_CLNDN.split('|')
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                jSZADiseaseID = []
                jDisease = []
                for k in range(0,len(kClinVarDiseases)):
                    kClinVarDisease_temp = kClinVarDiseases[k].replace(",", "").strip().lower()
                    matched = disease_lookup.get(kClinVarDisease_temp, None)
                    if matched:
                        did, dname, _ = matched
                        jVarSZAdiseaseIDs.append(did)
                        jVarDiseaseNames.append(dname)
                    else:
                        jSZADiseaseID.append('')  # 追加空字符串作为占位符
                        jDisease.append('') 
                jDisease = [x for x in jDisease if x != '']
                jSZADiseaseID = [x for x in jSZADiseaseID if x != '']
                if len(jDisease) == 1:
                    jDisease = ''.join(jDisease)
                    jSZADiseaseID = ''.join(jSZADiseaseID)
                    row = '\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases','Intergenic',jDisease]+iTemp)+'\n'
                    content_buffer.append(row)
                elif len(jDisease) > 1:
                    for m in range(len(jDisease)):
                        mjDisease = ''.join(jDisease[m])
                        mjSZADiseaseID = ''.join(jSZADiseaseID[m])
                        row = '\t'.join([jSZAvarID,jCLNID,jDatabaseType,mjSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases','Intergenic',mjDisease]+iTemp)+'\n'
                        content_buffer.append(row)

#%% Cell 8 Part B processing and generate outFile (saveFlag2) ==================================
# Part B
#1)MAF less than 0.05
#2)High Impact
#or
#2)
#ada_score>0.6
#rf_score >0.6
#REVEL>0.75
#Any 4 of the Splice_DS_score >0.5
#BayesDel_addAF_score >0.0692655
#BayesDel_noAF_score > -0.0570105
#newContentHeadings = ['SZAID','Database','Database_type','SZAdiseaseID','ClinVar_CLNSIG','ClinVar_CLNDN','SZAreportCategory','ClinVar_ID','selCSQ','Zygosity','Genotype']+geneBasedRefHeadings+headings;
#newContent = '\t'.join(newContentHeadings)+'\n'
COLUMN_INDEX = {
    'INFO': headings.index('INFO'),
    'CHROM': headings.index('#CHROM'),
    'REF': headings.index('REF'),
    'ALT': headings.index('ALT'),
    'SAMPLE': len(headings) - 1 if len(headings) > 9 else None
}

varCount = 0
for i in range(0,len(reportB)):
    iText = [s for s in reportB[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
    iTemp = reportB[i].replace('\n','').split('\t')
    iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
    ### ============================== Genotype processing part: if sample, get iGenotype and iZygo ===============
    # loop through each line of reportB
    if COLUMN_INDEX['SAMPLE'] is not None:
        sample_entry = reportB[i].split('\t')[COLUMN_INDEX['SAMPLE']]
        chrom = reportB[i].split('\t')[COLUMN_INDEX['CHROM']]
        ref = reportB[i].split('\t')[COLUMN_INDEX['REF']]
        alt = reportB[i].split('\t')[COLUMN_INDEX['ALT']]
        
        iGenotype, iZygo = gpf.parse_genotype(
            entry=sample_entry,
            ref=ref,
            alt=alt,
            chrom=chrom,
            check_male=check_male_flag
        )
    else:
        iGenotype = "Missing"
        iZygo = "Unknown"
    ### ============================== 
    varFlag = 0
    
    for j in range(0,len(iText)):
        jText = iText[j].split('|')
        pathogenicity = jText[col_map['ClinVar_CLNSIGCONF']].split('&')
        # Check MAX_AF threshold and known benign filter
        if (jText[col_map['MAX_AF']] == '' or float(jText[col_map['MAX_AF']]) < 0.05) and not any(s in pathogenicity for s in ['benign', 'likely_benign']):
            # Check pathogenicity impact conditions
            predicted_impact = (
                (jText[col_map['IMPACT']] == 'HIGH') or
                (jText[col_map['ada_score']] and float(jText[col_map['ada_score']]) > 0.6) or
                (jText[col_map['rf_score']] and float(jText[col_map['rf_score']]) > 0.6) or
                (jText[col_map['REVEL']] and float(jText[col_map['REVEL']]) > 0.75) or
                (jText[col_map['SpliceAI_pred_DS_AL']] and float(jText[col_map['SpliceAI_pred_DS_AL']]) > 0.5) or
                (jText[col_map['SpliceAI_pred_DS_DG']] and float(jText[col_map['SpliceAI_pred_DS_DG']]) > 0.5) or
                (jText[col_map['SpliceAI_pred_DS_DL']] and float(jText[col_map['SpliceAI_pred_DS_DL']]) > 0.5) or
                (jText[col_map['SpliceAI_pred_DS_AG']] and float(jText[col_map['SpliceAI_pred_DS_AG']]) > 0.5) or
                (jText[col_map['BayesDel_addAF_score']] and float(jText[col_map['BayesDel_addAF_score']]) > 0.0692655) or
                (jText[col_map['BayesDel_noAF_score']] and float(jText[col_map['BayesDel_noAF_score']]) > -0.0570105) or
                (jText[col_map['am_class']] == 'likely_pathogenic' and float(jText[col_map['am_pathogenicity']]) > 0.564)
            )
            if predicted_impact:
                jGenes = jText[3].split('&')
                if any([x for x in jGenes if x in geneBasedRef]):
                    varFlag = 1
                    jGeneList = [x for x in jGenes if x in geneBasedRef]
                    jInds = []
                    for k in range(0,len(jGeneList)):
                        kTemp = jGeneList[k]
                        jInds = jInds +[i for i, x in enumerate(geneBasedRef) if x == kTemp]
                    for k in range(0,len(jInds)):
                        scoreFlag = 2+(jText[col_map['IMPACT']]=='HIGH')*1
                        if ConfidenceLevel[jInds[k]] == 'High_confidence':
                            SZAscore = 10+scoreFlag
                        elif ConfidenceLevel[jInds[k]] == 'Moderate_confidence':
                            SZAscore = 5+scoreFlag
                        elif ConfidenceLevel[jInds[k]] == 'Low_confidence':
                            SZAscore = scoreFlag
                        else:
                            sys.exit('Unknown confidence score level for the report!')
                        jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                        #SZAscore = 2+('High_confidence' == ConfidenceLevel[jInds[k]])*5+(jText[colNames_CSQ.index('IMPACT')]=='HIGH')*1
                        try:
                            row = '\t'.join(
                            ['NovelTrans'+headings[9]+'_'+str(varCount+1),'NA','NA',
                                jSZADiseaseID,'NA','NA',str(SZAscore),'NA',iText[j],iZygo,
                                iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],
                                geneBasedRef[jInds[k]],Disease[jInds[k]]]+iTemp)+'\n'
                        except:
                            row = '\t'.join(
                            ['NovelTrans'+'_'+str(varCount+1),'NA','NA',
                                jSZADiseaseID,'NA','NA',str(SZAscore),'NA',iText[j],iZygo,
                                iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],
                                geneBasedRef[jInds[k]],Disease[jInds[k]]]+iTemp)+'\n'
                        content_buffer.append(row)
    if varFlag == 1:
        varCount = varCount+1

print("putative variants matched the GeneDB:",varCount)

# generate header line
header_line = '\t'.join(
    ['SZAID','Database','Database_type','SZAdiseaseID',
     'ClinVar_CLNSIG','ClinVar_CLNDN','SZAreportCategory',
     'ClinVar_ID','selCSQ','Zygosity','Genotype'] 
    + geneBasedRefHeadings 
    + headings
) + '\n'

# generete outFile
newContent = header_line + ''.join(content_buffer)
outfid = open(outFile,'w')
outfid.write(newContent)
outfid.close()
print("outFile for rare diseases information generated")


#%% Cell 9  get results file, spliting CSQ info and adding inheritance, disease info, remove duplicates ==============================
# Step A: read outFile
df_outFile = pd.read_csv(outFile, sep='\t',dtype={'ClinVar_ID': str})  # 相当于 final_report_sp
df_outFile['ClinVar_ID'] = df_outFile['ClinVar_ID'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
df_outFile['Database'] = df_outFile['Database'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)

# Step B: split 'selCSQ' => multiple columns
temp_df_outFile = df_outFile['selCSQ'].str.split('|', expand=True)
temp_df_outFile.columns = colNames_CSQ  # your existing column name list

# add columns that do not exist in df_outFile
new_columns = [col for col in temp_df_outFile.columns if col not in df_outFile.columns]
df_outFile = pd.concat([df_outFile, temp_df_outFile[new_columns]], axis=1)

# Step C: add Inheritances & DiseaseInfo
df_outFile['Inheritances'] = ""
df_outFile['DiseaseInfo']  = ""

for idx, row in df_outFile.iterrows():
    iSZAdiseaseID = row["SZAdiseaseID"]  # get the 4th column
    if 'SZA' not in str(iSZAdiseaseID):
        continue
    iGene = row["Genes"]  # get the 14th column
    iOMIM_num = DiseaseDB[DiseaseDB['SZAdiseaseID'] == iSZAdiseaseID]['DiseaseMIM'].to_list()[0]
    #print(iOMIM_num)
    if iOMIM_num:
        iheritance_list = OMIM_Inheritance_DB[
            (OMIM_Inheritance_DB['phenotypeMimNumber'] == iOMIM_num) &
            (OMIM_Inheritance_DB['approvedGeneSymbol'] == iGene)
        ]['inheritances'].to_list()
        if iheritance_list:
            iheritance = iheritance_list[0]
            #print(iheritance)
        else:
            iheritance = "No matched diseaseMIM in OMIM/OMIM not provided Inheritance"
    else:
        iheritance = ""

    iDiseaseInfo = (
        DiseaseDB[DiseaseDB['SZAdiseaseID'] == iSZAdiseaseID]['OMIM_Description'].to_list()[0]
        + '|' + "Reference:DiseaseName:"
        + DiseaseDB[DiseaseDB['SZAdiseaseID'] == iSZAdiseaseID]['DiseaseName'].to_list()[0]
        + '|' + "DiseaseSource:"
        + DiseaseDB[DiseaseDB['SZAdiseaseID'] == iSZAdiseaseID]['SourceName'].to_list()[0]
        + '|' + "Disease SourceID:"
        + DiseaseDB[DiseaseDB['SZAdiseaseID'] == iSZAdiseaseID]['SourceID'].to_list()[0]
        + '|' + "OMIM number:"
        + DiseaseDB[DiseaseDB['SZAdiseaseID'] == iSZAdiseaseID]['DiseaseMIM'].to_list()[0]
    )

    df_outFile.at[idx, 'Inheritances'] = iheritance
    df_outFile.at[idx, 'DiseaseInfo']  = iDiseaseInfo

# Step D: deduplication logic - keep the first line if (Target.group, Disease, Genes, SZAID, SZAreportCategory) are the same
df_outFile_dedup = df_outFile.drop_duplicates(
    subset=["Target.group", "Disease", "Genes", "SZAID", "SZAreportCategory"],
    keep="first"
).copy()

# Step E: add a column final_target_group at the beginning

df_outFile_dedup['final_target_group'] = df_outFile_dedup['Target.group'].apply(gpf.make_final_target_group)

# adjust the column order: first final_target_group, then the remaining columns
all_cols = list(df_outFile_dedup.columns)
all_cols.remove('final_target_group')
df_outFile_dedup = df_outFile_dedup[['final_target_group'] + all_cols]
cols_to_fix = ['Database', 'ClinVar_ID']  # Add any other affected columns
df_outFile_dedup[cols_to_fix] = df_outFile_dedup[cols_to_fix].astype(str)

# Step F: write the final nodup file

print("generating no dup file")

nodupfile = outFile.replace('.txt','_sp_Inheritance_4_nodup.txt')
df_outFile_dedup.to_csv(nodupfile, sep='\t', index=False)
#print(f"Generated no-dup file => {nodupfile}")

#%% Not RUN Cell 9 length file generation
# lengths = {
#     'reportA_ClinP_LP_var': len(reportA),
#     'reportB_Putative_var': len(reportB),
#     'reportB_Gene_matched':varCount,
#     'reportC_GWAS_var': len(reportC),
#     'reportD_Pharma_var': len(reportD),
#     'reportE_Trait_var': len(reportE),
#     'Gender': gender
# }
# # Create a DataFrame from the lengths dictionary
# df_lengths = pd.DataFrame(lengths.items(), columns=['Report Type', 'Length'])

# length_filename =outFile.replace('.txt','_length.txt') # Include fileName in the output filename

# df_lengths.to_csv(length_filename, index=False,sep='\t')

#%% Cell 10 GeneBe ACMG classification, ClinGen classification,GenCC and ClinGen gene-disease classification,ontology  ======================================
import genebe as gnb
genebe_df = pd.read_csv(nodupfile, sep="\t",header=0)

small_df = genebe_df.loc[:,["#CHROM","POS","REF","ALT"]]
small_df = small_df.rename(columns={"#CHROM":"chr","POS":"pos","REF":"ref","ALT":"alt"})
unique_small_df = small_df.drop_duplicates()

try:
    annotated_df = gnb.annotate(
        unique_small_df,
        genome='hg38',
        use_ensembl=False,
        use_refseq=True,
        flatten_consequences=True,
        output_format="dataframe"
    )
except Exception as e:
    print(f"GeneBe annotation failed: {str(e)}")
    # create an empty DataFrame to keep the structure
    annotated_df = pd.DataFrame(columns=[
        'chr', 'pos', 'ref', 'alt', 'gene_symbol', 
        'acmg_score', 'acmg_classification', 'acmg_criteria'
    ])
annotated_df = annotated_df.rename(columns={"chr":"#CHROM","pos":"POS","ref":"REF","alt":"ALT"})

required_columns = ["#CHROM","POS","REF","ALT","gene_symbol",
                   "acmg_score","acmg_classification",'acmg_criteria']
for col in required_columns:
    if col not in annotated_df.columns:
        annotated_df[col] = None  # add missing columns

small_annotate_all = annotated_df[required_columns].rename(columns={'gene_symbol': 'Genes'})
# merge automatically handles missing values
genebe_df = pd.merge(
    genebe_df, 
    small_annotate_all, 
    how="left", 
    on=["#CHROM","POS","REF","ALT","Genes"]
)
# the following code remains the same
genebe_df.to_csv(nodupfile, sep="\t", index=False)
print(genebe_df.shape)

genebe_df['ClinVar_ID'] = genebe_df['ClinVar_ID'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
genebe_df['Database'] = genebe_df['Database'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)



genebe_df["ClinVar_ID"] = genebe_df["Database"].str.split("&").str[0]

clingene = clingene[['#Variation', 'ClinVar Variation Id',
      'HGNC Gene Symbol', 'Disease', 'Mondo Id',
       'Mode of Inheritance', 'Assertion', 'Applied Evidence Codes (Met)',
       'Applied Evidence Codes (Not Met)', 'Summary of interpretation']]
clingene = clingene.rename(columns={"ClinVar Variation Id": "ClinVar_ID", "HGNC Gene Symbol":"Genes",
                                    "Assertion":"ClinGen_classification","Disease":"ClinGen_disease",
                                    "#Variation":"ClinGen_variant",
                                    "Applied Evidence Codes (Met)":"ClinGen_applied_evidence_codes"})
clingen_merge_df = pd.merge(genebe_df, clingene, how="left", on=["ClinVar_ID","Genes"])

print(clingen_merge_df.shape)

genecc_clingen_classification = pd.read_csv(geneBaseFile, sep="\t",header=0)
#genecc_clingen_classification.drop(columns='GenCC_classification_clingen', inplace=True)
clingen_merge_df = pd.merge(clingen_merge_df, genecc_clingen_classification, how="left",
                            on=['Gene.Disease.confidence.level', 'Target.group', 'Genes', 'Disease',
                                'SZAdiseaseID' ])

print(clingen_merge_df.shape)
ontology = ontology[['Genes',"SZAdiseaseID","category"]]
ontology = ontology.groupby(['Genes', 'SZAdiseaseID'])['category'].apply(lambda x: ', '.join(x)).reset_index()
clingen_merge_df_ontology = pd.merge(clingen_merge_df, ontology, on=['Genes','SZAdiseaseID'], how='left')
#clingen_merge_df.to_csv(nodupfile, sep="\t",index = None, header=True)
#col_relocation = ['final_target_group', 'SZAID', 'Database', 'Database_type', 'SZAdiseaseID', 'ClinVar_CLNSIG', 'acmg_classification', 'acmg_criteria', 'ClinGen_classification', 'ClinGen_applied_evidence_codes', 'ClinVar_CLNDN', 'SZAreportCategory', 'ClinVar_ID', 'selCSQ', 'Zygosity', 'Genotype', 'Gene.Disease.confidence.level', 'Target.group', 'Genes', 'Disease', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'REFSEQ_MATCH', 'REFSEQ_OFFSET', 'GIVEN_REF', 'USED_REF', 'BAM_EDIT', 'SOURCE', 'SIFT', 'PolyPhen', 'HGVS_OFFSET', 'AF', 'gnomADe_AF', 'gnomADe_AFR_AF', 'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF', 'gnomADe_SAS_AF', 'gnomADg_AF', 'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF', 'gnomADg_EAS_AF', 'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'ada_score', 'rf_score', 'REVEL', 'BayesDel_addAF_pred', 'BayesDel_addAF_rankscore', 'BayesDel_addAF_score', 'BayesDel_noAF_pred', 'BayesDel_noAF_rankscore', 'BayesDel_noAF_score', 'REVEL_rankscore', 'REVEL_score', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL', 'am_class', 'am_genome', 'am_pathogenicity', 'am_protein_variant', 'am_transcript_id', 'am_uniprot_id', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info', 'ClinVar', 'ClinVar_CLNHGVS', 'ClinVar_CLNSIGINCL', 'ClinVar_CLNVC', 'ClinVar_GENEINFO', 'ClinVar_CLNDISDB', 'ClinVar_CLNSIGCONF', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNDNINCL', 'Database_Type', 'Database_SZAID', 'Inheritances', 'DiseaseInfo', 'acmg_score', 'ClinGen_variant', 'ClinGen_disease', 'Mondo Id', 'Mode of Inheritance', 'Applied Evidence Codes (Not Met)', 'Summary of interpretation']
#clingen_merge_df = clingen_merge_df[col_relocation]
#clingen_merge_df = clingen_merge_df.fillna("")
clingen_merge_df_ontology['ClinVar_ID'] = clingen_merge_df_ontology['ClinVar_ID'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
clingen_merge_df_ontology['Database'] = clingen_merge_df_ontology['Database'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
clingen_merge_df_ontology.to_csv(nodupfile, sep="\t",index = None, header=True)

#%%cell 11 nodup4 file refinement  =======================================
# 1.HGVs naming converter
# Generate a nice HGVs format name on a new column "HGVS_Naming"
# format: NM_003126.4(SPTA1):c.6531-12C>T
# HGVSc + gene + HGVSp (if available)

# 2. HPO matching -- waiting
# generate a new column including HPO ontology names that matches the diseases

# 3.Known Benign variants filter
# based on CLIN_SIG and CLIVAR_CINCOF column , 
# filter out the known Benign variants for our conflicting variants and predicted variants

# 4.ClinVar Review star converter 
# For CLN Review Star, higher the star score, better confidence

# HGVS_Naming
clingen_merge_df_ontology['HGVS_Naming'] = clingen_merge_df_ontology.apply(gpf.generate_hgvs_naming, axis=1)

# ReviewStar
clingen_merge_df_ontology = gpf.add_review_star(clingen_merge_df_ontology)

# filter known benign variants
filtered_clingen_merge_df_ontology = gpf.filter_known_benign_variants(clingen_merge_df_ontology)
# process Genes column, only keep the first gene name
filtered_clingen_merge_df_ontology['Genes'] = filtered_clingen_merge_df_ontology['Genes'].str.split(';').str[0]

filtered_clingen_merge_df_ontology['ClinVar_ID'] = filtered_clingen_merge_df_ontology['ClinVar_ID'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
filtered_clingen_merge_df_ontology['Database'] = filtered_clingen_merge_df_ontology['Database'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)

# save results to new file
filtered_clingen_merge_df_ontology.to_csv(nodupfile, sep='\t', index=False)

print(filtered_clingen_merge_df_ontology.shape)
print(f"Refined nodup4 saved {nodupfile}")

# nodup files uploading
if Upload_To_Cloud == True:
	ssh = paramiko.SSHClient()
	ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	ssh.connect('ec2-16-171-6-238.eu-north-1.compute.amazonaws.com', username='ubuntu', password='AWuWqFRGHUhHV59!', port=22)
	scp = SCPClient(ssh.get_transport())
	remote_file_path = '/home/ubuntu/bckndsrv/RawDataOutputs(GenomicSZA)/' + nodupfile.split('/')[-1]
	scp.put(nodupfile, remote_file_path)
	print("File transferred successfully to AWS directory:", remote_file_path)

	# Close the SCP and SSH sessions
	scp.close()
	ssh.close()
#%%Cell 5 report C GWAS  newest update 20250210, speed up ---------------
start_time = time.time()
if run_GWAS == True:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportC)) + ' lines of GWAS variants' )
    # 1) Header prepare
    newLineHeadings = ['Trait','Genotype','Zygosity','risk_genotype','SZAvarID','GWASID','OR','Pval','PUBMED','Study accession ID','MAPPED_GENE',"CONTEXT"]
    newLine = '\t'.join(newLineHeadings) + '\n'
    gwas_outfile = outFile.replace('.txt','_GWAS.txt')
    with open(gwas_outfile, 'w') as f:
        f.write(newLine)
    # 2) to speed up, filter GWAS_db every time is too slow => pre-group by SZAID and store in a dictionary
    #    then directly O(1) get the corresponding line
    #    if GWAS_db is large, this will significantly speed up
    GWAS_dict = {}
    for szaid, group in GWAS_db.groupby('SZAID'):
        GWAS_dict[szaid] = group

    # 3) loop through each line of reportC
    for i in range(len(reportC)):
        line = reportC[i].rstrip('\n')
        fields = line.split('\t')  # only split once
        # iREF, iALT, genotype list
        iREF = fields[3].split(',')
        iALT = fields[4].split(',')
        iGenoTypeList = iREF + iALT
        # read sample genotype (assume in the 9th column)
        sample_geno_field = fields[9].split(':')[0]  # e.g. '0/1'
        # judge if missing
        if '.' not in sample_geno_field:
            iGenotypInd1 = int(sample_geno_field[0])
            iGenotypInd2 = int(sample_geno_field[2]) if len(sample_geno_field) > 2 else 0
            # try to construct iGenotype
            try:
                iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+ iGenoTypeList[iGenotypInd2]
            except:
                # original logic: fallback
                iGenotype = iREF[0] + '>' + iALT[0] + ':' + sample_geno_field
            # remove CSQ=... in INFO field (original logic)
            iTemp = fields[:]
            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
            # judge Zygosity
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
            # process GWAS annotation (only look at the first transcript)
            csq_part = [s for s in fields[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','')
            iText = csq_part.split(',')
            jText = iText[0].split('|')
            # judge if Database_Type contains GWAS_var
            if jText[col_map['Database_Type']] != '':
                db_type_vals = jText[col_map['Database_Type']].split('&')
                if 'GWAS_var' in db_type_vals:
                    # find the index of GWAS_var
                    jGWASind = db_type_vals.index('GWAS_var')
                    jGWASID = jText[col_map['Database']].split('&')[jGWASind]
                    jSZAvarID = jText[col_map['Database_SZAID']].split('&')[jGWASind]
                    # 4) directly get GWAS_db_filtered from the dictionary (for speed up)
                    GWAS_db_filtered = GWAS_dict.get(jSZAvarID, None)
                    if GWAS_db_filtered is not None:
                        # loop through each record
                        for k, row in GWAS_db_filtered.iterrows():
                            jOR = row['OR.or.BETA']
                            jPval = row['P.VALUE']
                            jrisk_genotype = row['risk_genotype']
                            jPUBMED = row['PUBMEDID']
                            jTrait = row['DISEASE.TRAIT']
                            jMappedGene = row['MAPPED_GENE']
                            jVarType = row['CONTEXT']
                            jStudy_accession = row['STUDY.ACCESSION']
                            # only keep p <= 5e-8 and iGenotype == risk_genotype (keep original logic)
                            if iGenotype == jrisk_genotype and jPval <= 5e-8:
                                newLine_temp = '\t'.join([
                                    jTrait,
                                    iGenotype,
                                    iZygo,
                                    jrisk_genotype,
                                    jSZAvarID,
                                    jGWASID,
                                    str(jOR),
                                    str(jPval),
                                    str(jPUBMED),
                                    str(jStudy_accession),
                                    jMappedGene,
                                    jVarType
                                ]) + '\n'
                                # write to the same output file
                                with open(gwas_outfile, 'a') as f:
                                    f.write(newLine_temp)
        # print progress every 1000 lines
        if i % 10000 == 0:
            print(str(i) + ' lines processed!')

    end_time = time.time()
    print("Total processing time: {:.2f} seconds".format(end_time - start_time))
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

    if Upload_To_Cloud == True:
        # Initialize SSH client for file transfer
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect('ec2-16-171-6-238.eu-north-1.compute.amazonaws.com', username='ubuntu', password='AWuWqFRGHUhHV59!', port=22)
        scp = SCPClient(ssh.get_transport())

        # Define the local and remote paths
        local_file_path = outFile.replace('.txt', '_GWAS.txt')
        remote_file_path = '/home/ubuntu/bckndsrv/RawDataOutputs(GenomicSZA)/' + local_file_path.split('/')[-1]

        # Transfer the file
        scp.put(local_file_path, remote_file_path)
        print("File transferred successfully to AWS directory.")

        # Close SCP and SSH sessions
        scp.close()
        ssh.close()

#%% Cell 6 report D for pharmaco-annotation newest 20250210
# adding rsID after SZAvarID 20250110
start_time = time.time()
if run_Pharmaco == True:
    now = datetime.now()
    newLineHeadings = [
        'Drug.s.', 'Genotype', 'Zygosity', 'pharma_genotype', 'SZAvarID', 'rsID',
        'Annotation.Text', 'Score', 'Phenotype.Category', 'Phenotype.s.',
        'Gene', 'Level.of.Evidence', 'Level.Modifiers', 'URL'
    ]
    newLine = '\t'.join(newLineHeadings) + '\n'
    out_pharma = outFile.replace('.txt','_pharmaco.txt')
    with open(out_pharma, 'w') as f:
        f.write(newLine)

    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportD)) + ' lines of pharmaco variants' )

    # === 1) to speed up, pre-group Pharma_db by SZAID and store in a dictionary ===
    Pharma_dict = {}
    for szaid, group in Pharma_db.groupby('SZAID'):
        Pharma_dict[szaid] = group

    for i in range(len(reportD)):
        line = reportD[i].rstrip('\n')
        fields = line.split('\t')
        iREF = fields[3].split(',')
        iALT = fields[4].split(',')
        iGenoTypeList = iREF + iALT

        # judge if genotype is missing
        sample_geno_field = fields[9].split(':')[0]  # e.g. '0/1'
        if '.' not in sample_geno_field:
            iGenotypInd1 = int(sample_geno_field[0])
            iGenotypInd2 = int(sample_geno_field[2]) if len(sample_geno_field) > 2 else 0

            # construct genotype, main logic + fallback
            try:
                iGenotype = iGenoTypeList[0] + '/' + iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] + '/' + iGenoTypeList[iGenotypInd2]
            except:
                iGenotype = iREF[0] + '>' + iALT[0] + ':' + sample_geno_field

            # remove CSQ in INFO field
            fields_copy = fields[:]
            fields_copy[headings.index('INFO')] = fields_copy[headings.index('INFO')].split(';CSQ=')[0]

            # judge Zygosity
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
            # extract annotation of the first transcript
            info_part = fields[7]
            csq_candidates = [s for s in info_part.split(';') if 'CSQ=' in s]
            if csq_candidates:
                iText = csq_candidates[0].replace('CSQ=', '').split(',')
                jText = iText[0].split('|')
                if jText[col_map['Database_Type']] != '':
                    db_type_list = jText[col_map['Database_Type']].split('&')
                    if 'Pharma_var' in db_type_list:
                        jPharmaind = db_type_list.index('Pharma_var')
                        jPharmaID   = jText[col_map['Database']].split('&')[jPharmaind]
                        jSZAvarID   = jText[col_map['Database_SZAID']].split('&')[jPharmaind]

                        # === 2) use the pre-built Pharma_dict to get the sub-table ===
                        Pharma_db_filtered = Pharma_dict.get(jSZAvarID, None)
                        if Pharma_db_filtered is not None:
                            # old logic: only find "first" match iGenotype
                            # => iGenotype in jpharma_genotype_list
                            jpharma_genotype_list = Pharma_db_filtered['pharma_genotype'].tolist()
                            if iGenotype in jpharma_genotype_list:
                                genotype_ind = jpharma_genotype_list.index(iGenotype)
                                # get the corresponding row
                                match_row = Pharma_db_filtered.iloc[genotype_ind]
                                jpharma_genotype = match_row['pharma_genotype']
                                jDrug            = match_row['Drug.s.']
                                jrsID            = match_row['ID']
                                jAnnotation      = match_row['Annotation.Text']
                                jScore           = match_row['Score']
                                jPheno_category  = match_row['Phenotype.Category']
                                jPheno           = match_row['Phenotype.s.']
                                jGene            = match_row['Gene']
                                jevidence        = match_row['Level.of.Evidence']
                                jlevel_modifier  = match_row['Level.Modifiers']
                                jURL             = match_row['URL']

                                newLine_temp = '\t'.join([
                                    str(jDrug),
                                    iGenotype,
                                    iZygo,
                                    str(jpharma_genotype),
                                    jSZAvarID,
                                    str(jrsID),
                                    str(jAnnotation),
                                    str(jScore),
                                    str(jPheno_category),
                                    str(jPheno),
                                    str(jGene),
                                    str(jevidence),
                                    str(jlevel_modifier),
                                    str(jURL)
                                ]) + '\n'
                                # append to the same output file
                                with open(out_pharma, 'a') as f:
                                    f.write(newLine_temp)

        # print progress
        if i % 1000 == 0:
            print(str(i) + " lines processed")

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    end_time = time.time()
    print("Total processing time: {:.2f} seconds".format(end_time - start_time))
    if Upload_To_Cloud == True:
    # Initialize SSH client for file transfer
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect('ec2-16-171-6-238.eu-north-1.compute.amazonaws.com', username='ubuntu', password='AWuWqFRGHUhHV59!', port=22)
        scp = SCPClient(ssh.get_transport())
        # Define the remote path
        local_file_path = outFile.replace('.txt', '_pharmaco.txt')
        remote_file_path = '/home/ubuntu/bckndsrv/RawDataOutputs(GenomicSZA)/' + local_file_path.split('/')[-1]
        # Transfer the file
        scp.put(local_file_path, remote_file_path)
        print("File transferred successfully to AWS directory.")

        # Close SCP and SSH sessions
        scp.close()
        ssh.close()

#%% Cell 7 Traits report ====================================
## note: one SNP (rsID) could match to more than one trait phenotypes. such as rs1815739. In the result file we only output one db matching but we do the multiple match in the updated_traits.csv
if run_traits == True:
    from datetime import datetime
    now = datetime.now()
    newLineHeadings = ['Traits name', 'category', 'genes', 'variants', 'Description','Patient genotypes','Genotype Description','Zygosity'];
    newLine = '\t'.join(newLineHeadings)+'\n'
    with open(outFile.replace('.txt','_traits.txt'),'w') as f:
        f.write(newLine)
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportE)) + ' lines of traits' )
    #this irs_id_set is used to record the trait variants this individual have, and later write the one don't have to a "wild type"
    irs_id_set = {s.split('|')[col_map['Existing_variation']].split('&')[0] for s in reportE}
    for i in range(0,len(reportE)):
        current_time = now.strftime("%H:%M:%S")
        ##iText means all transcripts
        iText = [s for s in reportE[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
        irs_id = iText[0].split("|")[col_map['Existing_variation']].split('&')[0]
        iREF = reportE[i].split('\t')[3].split(',')
        iALT = reportE[i].split('\t')[4].split(',')
        iGenoTypeList = iREF+iALT
        if "." not in reportE[i].split('\t')[9].split(':')[0]:
            iGenotypInd1 = int(reportE[i].split('\t')[9][0])
            iGenotypInd2 = int(reportE[i].split('\t')[9][2])
            iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            #try:
            #    iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            #except:
            #    iGenotype = iREF[0]+'>'+iALT[0]+':'+iGenoType+reportA[i].split('\t')[9][0]+reportA[i].split('\t')[9][1]+reportA[i].split('\t')[9][2]
            iTemp = reportE[i].replace('\n','').split('\t')
            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
            #sys.exit('!')
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
            ##get the Pharma and sza id from sub transcripts, named jText
            ## here we seems not care about genes, transcripts, only focuse on variant so I deleted for loop
            iTraits= Trait_db[Trait_db['variants']== irs_id]['Traits name'].tolist()[0]
            iCategory = Trait_db[Trait_db['variants']== irs_id]['category'].tolist()[0]
            iGenes= Trait_db[Trait_db['variants']== irs_id]['genes'].tolist()[0]
            iDescription = Trait_db[Trait_db['variants']== irs_id]['Description'].tolist()[0].strip()
            iGenoType_description =Trait_db[Trait_db['variants']== irs_id]['Genotype_Description'].tolist()[0].strip()
            newLine_temp ='\t'.join([iTraits,iCategory,iGenes,irs_id,iDescription,iGenotype,iGenoType_description,iZygo])+'\n'
            with open(outFile.replace('.txt','_traits.txt'),'a') as f:
                f.write(newLine_temp)
    # Iterate through Trait_db and write entries not in reportE
    for index, row in Trait_db.iterrows():
        if row['variants'] not in irs_id_set:
            print(row['variants'] + "This patient carries wildtype genotype (keep as homo_ref)")
            newLine_temp = '\t'.join([row['Traits name'], row['category'], row['genes'], row['variants'], row['Description'], row['REF']+'/'+row['REF'], row['Genotype_Description'], 'homo_ref']) + '\n'
            with open(outFile.replace('.txt', '_traits.txt'), 'a') as f:
                f.write(newLine_temp)
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time) 

#%% APOE haploytype based on generated results file
    
rows = []
trait_outfile = outFile.replace('.txt','_traits.txt')

#  open the original file, read line by line
with open(trait_outfile, 'r', encoding='utf-8') as f_in:
    lines = f_in.readlines()

rs429358_gt = None
rs7412_gt = None

for line in lines:
    parts = line.strip().split('\t')
    if len(parts) < 6:
        continue  # skip incomplete lines
    
    # parts[3] => 'rs429358' or 'rs7412'
    variant = parts[3]
    # parts[5] => 'Patient genotypes' (e.g. "T/T", "C/C")
    genotype_str = parts[5]
    if '>' in genotype_str:
        genotype = genotype_str.split('>')[1]  # get the part after '>'
    else:
        genotype = genotype_str  # no '>', use the original string directly

    genotype = genotype.replace('/', '')  
    if variant == 'rs429358':
        rs429358_gt = genotype
        print(parts[5])
        print("rs429358_gt",rs429358_gt)
    elif variant == 'rs7412':
        rs7412_gt = genotype
        print("rs7412",rs7412_gt)

# 3) if both are got, then calculate Isoform & Risk
if rs429358_gt and rs7412_gt:
    print("get APOE genotype from original trait file")
    isoform, risk = gpf.get_apoe_isoform_and_risk(rs429358_gt, rs7412_gt)
    newLine_temp_parts = [
        "Alzheimer's disease risk",                     # Traits name
        "Neurogenic and Cognitive functions",          # category
        "APOE",                                        # genes
        "rs429358+rs7412",                               # variants
        "Combined APOE genotype from rs429358+rs7412", # Description
        isoform,                  # Patient genotypes
        risk,                                       # Genotype Description
        ""                                           # Zygosity
    ]
    newLine_temp = '\t'.join(newLine_temp_parts) + '\n'

    with open(trait_outfile, 'a', encoding='utf-8') as f_out:
        f_out.write(newLine_temp)
    print(f"Done. Check: {trait_outfile}")
else:
    isoform, risk = ("Unknown", "Unknown")
    print("did not find two genptypes of rs429358 and rs7412")

    

if Upload_To_Cloud == True:
# Initialize SSH client for file transfer
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('ec2-16-171-6-238.eu-north-1.compute.amazonaws.com', username='ubuntu', password='AWuWqFRGHUhHV59!', port=22)
    scp = SCPClient(ssh.get_transport())

    # Define the local and remote paths
    local_file_path = outFile.replace('.txt', '_traits.txt')
    remote_file_path = '/home/ubuntu/bckndsrv/RawDataOutputs(GenomicSZA)/' + local_file_path.split('/')[-1]
    # Transfer the file
    scp.put(local_file_path, remote_file_path)
    print("File transferred successfully to AWS directory.")
    # Close SCP and SSH sessions
    scp.close()
    ssh.close()


all_end_time = time.time()
print("All Process total processing time: {:.2f} seconds".format(all_end_time - all_start_time))

# 预处理疾病数据库
disease_map = {
    name.strip().lower().replace(',', ''): (did, name) 
    for did, name in zip(DiseaseID_DSDB, DiseaseName_DSDB)
}

def match_disease(clinvar_disease: str) -> tuple:
    key = clinvar_disease.strip().lower().replace(',', '')
    return disease_map.get(key, ('', clinvar_disease))