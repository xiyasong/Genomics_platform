#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:41:09 2024

@author: xinmengliao
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:41:09 2023

@author: xiyas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 17:40:59 2023

@author: xiyas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 11:31:41 2022

@author: xiyas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 01:33:22 2022

@author: xiyas
"""

# GenomeReportingPipeline version 2 on SZA sercer
# June 2023: Adding the male/female distinguishment for "Hemizygous"
# Jan 2024: vep gz file compatible

#%%Cell 1 Define the path for parameters and databases
import sys
import pandas as pd
import numpy as np
import os
import gzip
import re

# 1. choose which enviroment running local uppmax, SZA
path = 'server'

# 2. choose which cohort belongs to, Turkish or Swedish
cohort = 'Turkish'
print("Turkish vcf files")
# 3. choose whether need to run the GWAS and Pharmaco findings
run_GWAS = True
run_Pharmaco = True
run_traits = True
run_simplify = False
Upload_To_Cloud = False

# Check the input value and execute code accordingly
if path == 'local':
    fileName = sys.argv[1]
    #fileName = '/Users/xinmengliao/Documents/Project/20231227_Newborn/Scripts/P0064_821_vep_annotated.vcf.gz'
    #fileName = '/Volumes/SSamantha/Newborn_843_v1/vep/P0064_821_vep_annotated.vcf.gz'
    #output_directory = sys.argv[2]
    outFile = sys.argv[2]
    #outFile = '/Users/xinmengliao/Documents/Project/20231227_Newborn/Scripts/P0064_821_vep_annotated.txt'

    # # define name for LOFTEE annotated file
    # if loftee_plugin == True:
    #     print("LOFTEE plugin results will be included in the final report.")
    #     loftee = sys.argv[3]    
    # #loftee = '/Users/xinmengliao/Documents/Project/20231227_Newborn/Scripts/P0064_1_loftee_vep_annotated.vcf.gz'

    # ##define the name for the output file
    geneBaseFile = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Genedb1213.txt'
    diseaesDBfile = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/diseaseDB_0605.txt'
    OMIM_inheritance_DBfile = '/Users/xinmengliao/Documents/Project/20231015_SZA/database-file-v2/pheno_OMIM_all.txt'
    ontology_file = "/Users/xinmengliao/Documents/Project/20231227_Newborn/Ontology/genedb.ontology.all0307.csv"

    #Read disease and inheritance DB
    DiseaseDB = pd.read_csv(diseaesDBfile, sep="\t",encoding="ISO-8859-1")
    DiseaseDB= DiseaseDB.replace(np.nan,"No info")
    OMIM_Inheritance_DB = pd.read_csv(OMIM_inheritance_DBfile,sep="\t",dtype={"phenotypeMimNumber": str})
    OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")
    HGMD_file = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Python pipeline/HGMD_DM20240917.txt'
    HGMD_DM = pd.read_csv(HGMD_file, sep="\t")
    HGMD_chrom_dict = {}

    for index, row in HGMD_DM.iterrows():
        variant_info = row['variant_info']
        chrom, pos_ref_alt = variant_info.split('_', 1)
    
        if chrom not in HGMD_chrom_dict:
            HGMD_chrom_dict[chrom] = set()
    
        HGMD_chrom_dict[chrom].add(pos_ref_alt)

    #Read omim, hgnc, decipher files
    hgnc_file = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Python pipeline/hgnc.txt'
    omim_file = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Python pipeline/omim.txt'
    decipher_file = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Python pipeline/decipher.txt'
    clingen_file = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Clingen/Clingen-variants-2024-12-09.txt'

    print("Set path for local")

    if run_GWAS == True:
        GWAS_dbfile = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Python pipeline/Merged_GWAS_vcf0627.txt'
        GWAS_db =pd.read_csv(GWAS_dbfile, sep="\t")
        GWAS_db= GWAS_db.replace(np.nan,"No info")
        print("GWAS matching will run")

    if run_Pharmaco == True:
        #Pharma_dbfile = '/Users/xinmengliao/Documents/Project/20231015_SZA/database-file-v2/Merged_Pharma_vcf.txt'
        Pharma_dbfile = '/Users/xinmengliao/Documents/Project/20231227_Newborn/PharmGKB/AllPGx_annotation1218.txt'
        Pharma_db =pd.read_csv(Pharma_dbfile, sep="\t")
        #Pharma_db= Pharma_db.replace(np.nan,"No info")
        # Read the haplotype variants data
        HaplotypeFile='/Users/xinmengliao/Documents/Project/20231227_Newborn/PharmGKB/All_Haplotype_var1218.txt'
        with open(HaplotypeFile, 'r') as file:
            haplotype_var = file.read().strip().splitlines()
        HaplotypeID='/Users/xinmengliao/Documents/Project/20231227_Newborn/PharmGKB/All_Haplotype_rsID1218.txt'
        with open(HaplotypeID, 'r') as file:
            haplotype_rsID = file.read().strip().splitlines()
        print("Pharmaco matching will run")

    if run_traits == True:
        Trait_dbfile = '/Users/xinmengliao/Documents/Project/20231227_Newborn/database_files/Python pipeline/Reports_genome_databases_traits_merged.txt'
        Trait_db =pd.read_csv(Trait_dbfile,sep="\t",encoding = "ISO-8859-1")
        Trait_db = Trait_db.dropna(axis=1, how='all')
        Trait_db= Trait_db.replace(np.nan,"No info")
        trait_list=Trait_db['variants'].to_list()
        print("Traits matching will run")
    
    if run_simplify:
        print("Simplify result file witll run")


elif path == 'server':
    fileName = sys.argv[1]
    outFile = sys.argv[2]

    # ##define the name for the output file
    geneBaseFile = '/mnt/storage_pool/Genomics/Genome/database-files/Genedb1213.txt'
    diseaesDBfile = '/mnt/storage_pool/Genomics/Genome/database-files/diseaseDB_0605.txt'
    OMIM_inheritance_DBfile = '/mnt/storage_pool/Genomics/Genome/database-files/pheno_OMIM_all.txt'
    ontology_file = '/mnt/storage_pool/Genomics/Genome/database-files/genedb.ontology.all0307.csv'

    #Read disease and inheritance DB
    DiseaseDB = pd.read_csv(diseaesDBfile, sep="\t",encoding="ISO-8859-1")
    DiseaseDB= DiseaseDB.replace(np.nan,"No info")
    OMIM_Inheritance_DB = pd.read_csv(OMIM_inheritance_DBfile,sep="\t",dtype={"phenotypeMimNumber": str})
    OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")
    HGMD_file = '/mnt/storage_pool/Genomics/Genome/database-files/HGMD_DM20240917.txt'
    HGMD_DM = pd.read_csv(HGMD_file, sep="\t")
    HGMD_chrom_dict = {}

    for index, row in HGMD_DM.iterrows():
        variant_info = row['variant_info']
        chrom, pos_ref_alt = variant_info.split('_', 1)
    
        if chrom not in HGMD_chrom_dict:
            HGMD_chrom_dict[chrom] = set()
    
        HGMD_chrom_dict[chrom].add(pos_ref_alt)

    #Read omim, hgnc, decipher files
    hgnc_file = '/mnt/storage_pool/Genomics/Genome/database-files/hgnc.txt'
    omim_file = '/mnt/storage_pool/Genomics/Genome/database-files/omim.txt'
    decipher_file = '/mnt/storage_pool/Genomics/Genome/database-files/decipher.txt'
    clingen_file = '/mnt/storage_pool/Genomics/Genome/database-files/Clingen-variants-2024-12-09.txt'

    print("Set path for SZA")

    if run_GWAS == True:
        GWAS_dbfile = '/mnt/storage_pool/Genomics/Genome/database-files/Merged_GWAS_vcf_2024.txt'
        GWAS_db =pd.read_csv(GWAS_dbfile, sep="\t")
        GWAS_db= GWAS_db.replace(np.nan,"No info")
        print("GWAS matching will run")

    if run_Pharmaco == True:
        Pharma_dbfile = '/mnt/storage_pool/Genomics/Genome/database-files/AllPGx_annotation1218.txt'
        Pharma_db =pd.read_csv(Pharma_dbfile, sep="\t")
        #Pharma_db= Pharma_db.replace(np.nan,"No info")
        # Read the haplotype variants data
        HaplotypeFile='/mnt/storage_pool/Genomics/Genome/database-files/All_Haplotype_var1218.txt'
        with open(HaplotypeFile, 'r') as file:
            haplotype_var = file.read().strip().splitlines()
        HaplotypeID='/mnt/storage_pool/Genomics/Genome/database-files/All_Haplotype_rsID1218.txt'
        with open(HaplotypeID, 'r') as file:
            haplotype_rsID = file.read().strip().splitlines()
        print("Pharmaco matching will run")

    if run_traits == True:
        Trait_dbfile = '/mnt/storage_pool/Genomics/Genome/database-files/Reports_genome_databases_traits_merged.txt'
        Trait_db =pd.read_csv(Trait_dbfile,sep="\t",encoding = "ISO-8859-1")
        Trait_db = Trait_db.dropna(axis=1, how='all')
        Trait_db= Trait_db.replace(np.nan,"No info")
        trait_list=Trait_db['variants'].to_list()
        print("Traits matching will run")



##NEW update 2.07: I add a new annotation, Clinvar review status for varisnts. so the CSQ from 88 to 89 columns
# Need to add new ClinVar disease index col if using the new ClinVar vcf.gz file (20240603)

# When using --refseq in vep annotation, more columns exist in the vep_files. 

colNames_CSQ =['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp',
'cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation',
'DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','REFSEQ_MATCH',
'REFSEQ_OFFSET','GIVEN_REF','USED_REF','BAM_EDIT','SOURCE','SIFT','PolyPhen','HGVS_OFFSET','AF',
'gnomADe_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_ASJ_AF','gnomADe_EAS_AF','gnomADe_FIN_AF',
'gnomADe_MID_AF','gnomADe_NFE_AF','gnomADe_REMAINING_AF','gnomADe_SAS_AF','gnomADg_AF','gnomADg_AFR_AF',
'gnomADg_AMI_AF','gnomADg_AMR_AF','gnomADg_ASJ_AF','gnomADg_EAS_AF','gnomADg_FIN_AF','gnomADg_MID_AF',
'gnomADg_NFE_AF','gnomADg_REMAINING_AF','gnomADg_SAS_AF','MAX_AF','MAX_AF_POPS','CLIN_SIG','SOMATIC','PHENO',
'ada_score','rf_score','REVEL','BayesDel_addAF_pred','BayesDel_addAF_rankscore','BayesDel_addAF_score',
'BayesDel_noAF_pred','BayesDel_noAF_rankscore','BayesDel_noAF_score','REVEL_rankscore','REVEL_score','SpliceAI_pred_DP_AG',
'SpliceAI_pred_DP_AL','SpliceAI_pred_DP_DG','SpliceAI_pred_DP_DL','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG',
'SpliceAI_pred_DS_DL','SpliceAI_pred_SYMBOL','am_class','am_genome','am_pathogenicity','am_protein_variant','am_transcript_id',
'am_uniprot_id','LoF','LoF_filter','LoF_flags','LoF_info','ClinVar','ClinVar_ID','ClinVar_CLNSIG','ClinVar_CLNDN','ClinVar_CLNHGVS',
'ClinVar_CLNSIGINCL','ClinVar_CLNVC','ClinVar_GENEINFO','ClinVar_CLNDISDB','ClinVar_CLNSIGCONF','ClinVar_CLNREVSTAT',
'ClinVar_CLNDNINCL','Database','Database_Type','Database_SZAID']




#%%Cell 2 start the analysis program for the vep annotated vcf files, generated report A B C and D || E for the traits
#file_count = 0
#for name in os.listdir("/proj/snic2020-16-69/nobackup/WGS_SZA/vep_annotated_vcf"):
#    print(name)
#    file_count +=1
#    print(file_count)
#    fileName = os.path.join("/proj/snic2020-16-69/nobackup/WGS_SZA/vep_annotated_vcf", name)
    #fileName = sys.argv[1]
    #geneBaseFile = sys.argv[2]
    #diseaesDBfile= sys.argv[3]
    #OMIM_inheritance_DBfile = sys.argv[4]
#    outFile = name + 'outfile.txt'
    #fileName = '/Users/xiyas/Report_genome/Sample_data_SZA/anno_P001_106_new_0830.vcf'
    #outFile = '/Users/xiyas/Report_genome/Sample_data_SZA/finalReport_P001_106_python.txt'
# read VEP output


check_male_flag = False
file = gzip.open(fileName,'rt')
tLine = file.readline()
#tLine = tLine.decode("ISO-8859-1")
i = 0
reportA ,reportB,reportC,reportD,reportE = [], [], [], [], []
while tLine.endswith('\n'):
    i = i+1
    iContent = tLine.replace('\n','').split('\t')
    #####check it's male or female, to get the correct Zygosity for the chrX variants
    #check_male_flag = 0>1
    ##get the content from VCF annotation header
    if '#' in iContent[0]:
        if 'ID=CSQ' in tLine:
            annoText = iContent[0].split('Format: ')
            colNames_CSQ = annoText[1].replace('">','')
            colNames_CSQ = colNames_CSQ.split('|')
        elif iContent[0] == '#CHROM':
            headings = iContent
    #start processing
    else:
        iText = [s for s in iContent[headings.index('INFO')].split(';') if 'CSQ=' in s]
        iText = iText[0].replace('CSQ=','').split(',')
        saveFlag1, saveFlag2,saveFlag3,saveFlag4, saveFlag5,saveFlag6 = False, False, False, False, False, False

        ##this is going to analyse all transcript result and any one satisfied the criteria, the whole variant goes into report A/B.

        ## 1.26 2024: In the conference I checked rs6025, a famous PGx genes for F5 gene. However it only appears on nodup4files, as a clinVar genes,
        ## not in PGx reports. so I realized this elif here should be all changed to "if", because elif ignores other cases if the first one satisfied the criteria
        for j in range(0,len(iText)):
            jText = iText[j].split('|')

            ichr = iContent[headings.index('#CHROM')].split('_')[0]
            ipos = iContent[headings.index('POS')]
            iref = iContent[headings.index('REF')]
            ialt = iContent[headings.index('ALT')]
            ivariation3 = f"{ipos}_{iref}_{ialt}"
            ivariation4 = f"{ichr}_{ipos}_{iref}_{ialt}"

            is_clinvar = jText[colNames_CSQ.index('Database_SZAID')] != '' and 'ClinP_LP_var' in jText[colNames_CSQ.index('Database_Type')].split('&')

            if is_clinvar:
                saveFlag1 = 'ClinP_LP_var'

            if ichr in HGMD_chrom_dict and ivariation3 in HGMD_chrom_dict[ichr] and is_clinvar == False :
                saveFlag6 = 'HGMD_var'
                
            #if not saveFlag1 and ((jText[colNames_CSQ.index('MAX_AF')] == '') or (float(jText[colNames_CSQ.index('MAX_AF')]) < 0.05)):
            #2024.0311 update: gnomad Version 4.0 allele frequencies by VEP annotation
            
            ## step 1 : get the MAX_AF --------------------------------------
            #start_index = colNames_CSQ.index('gnomAD4_AF')
            #end_index = colNames_CSQ.index('gnomAD4_AF_sas') + 1 
            if jText[colNames_CSQ.index('MAX_AF')] == '' : 
                jText[colNames_CSQ.index('MAX_AF')] = 0
            jText[colNames_CSQ.index('MAX_AF')] = float(jText[colNames_CSQ.index('MAX_AF')])
            #gnomad_values = [float(val) if val != '' else 0 for val in jText[start_index:end_index]]
            #max_gnomad_value = max(gnomad_values)
            #max_af_value = jText[colNames_CSQ.index('MAX_AF')]
            #max_value = max(float(max_gnomad_value), float(max_af_value))
            #jText[colNames_CSQ.index('MAX_AF')] = max_value
            #print(max_af_value)
            #print(max_gnomad_value)
            ## step 2: use the new MAX_AF from the comparison
            # ------------------------------------------------------------
            if not saveFlag1 and not saveFlag6 and jText[colNames_CSQ.index('MAX_AF')] < 0.05:
                ##Update: to clear the logics here need to add brackets
                if (jText[colNames_CSQ.index('IMPACT')] == 'HIGH' or (
                        jText[colNames_CSQ.index('ada_score')] != '' and float(jText[colNames_CSQ.index('ada_score')])>0.6) or (
                        jText[colNames_CSQ.index('rf_score')] != '' and float(jText[colNames_CSQ.index('rf_score')])>0.6) or (
                        jText[colNames_CSQ.index('REVEL')]!= '' and float(jText[colNames_CSQ.index('REVEL')])>0.75) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')])>0.5) or (
                        jText[colNames_CSQ.index('BayesDel_addAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_addAF_score')])>0.0692655) or (
                        jText[colNames_CSQ.index('BayesDel_noAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_noAF_score')])>-0.0570105) or (
                        jText[colNames_CSQ.index('am_class')]== 'likely_pathogenic' and float(jText[colNames_CSQ.index('am_pathogenicity')])> 0.564) or (
                        jText[colNames_CSQ.index('LoF')]== 'HC')):
                    saveFlag2 ="Putative_var"
                    #MAX_AF to gnomad global AF, and new filters need to be applied here too
            
            if jText[colNames_CSQ.index('Database_SZAID')] != '' and 'GWAS_var' in jText[colNames_CSQ.index('Database_Type')].split('&'):
                saveFlag3 = "GWAS_var"
            if jText[colNames_CSQ.index('Database_SZAID')] != '' and 'Pharma_var' in jText[colNames_CSQ.index('Database_Type')].split('&'):
                saveFlag4 = "Pharma_var"

            if run_Pharmaco == True:  
            #Check if Haplotype location is in the file 
                if ivariation4 in haplotype_var:
                    saveFlag4 = "Pharma_var"    

                #Check if Haplotype rsID is in the file 
                rsID = jText[colNames_CSQ.index('Existing_variation')]
                rsID = rsID.split('&')
                if True in [part in haplotype_rsID for part in rsID]:
                    saveFlag4 = "Pharma_var" 
            
            if run_traits == True:
                if jText[colNames_CSQ.index('Existing_variation')] in trait_list:
                    saveFlag5 ="Trait_var"
            
            # Check gender 
            if 'chrY' in iContent[headings.index('#CHROM')] :
                check_male_flag = True

        if saveFlag1:
            #if len(jText[colNames_CSQ.index('Database_Type')].split('&'))>1:
            #    sys.exit('!')
            reportA.append(tLine)
        if saveFlag2:
            reportB.append(tLine)
        if saveFlag3:
            reportC.append(tLine)
        if saveFlag4:
            reportD.append(tLine)
        if saveFlag5:
            reportE.append(tLine)
        if saveFlag6:
            tLine = "HGMD_" + tLine
            reportA.append(tLine)


#    if "chrY" in tLine:
#            #print("Sex:Male")
#        check_male_flag ="Male"
    
    #    if "27771923" in tLine:
#        sys.exit("!")

    if i%100000 == 0:
        print(str(i)+' lines processed!')
    tLine = file.readline()
    
variant_count = i
gender = "Male" if check_male_flag else "Female"

if gender == "Male" :
    print("Male")
else:
    print("Female")

file.close()
print('Cell 2 VEP annotated File processing done! Now start to map GeneDB and DiseaseDB')



#%%Cell 3 read gene db and disease db

print("running cell 3,reading gene db and disease db")
# read gene based database
ConfidenceLevel = []
TargetGroup = []
geneBasedRef = []
Disease = []
SZAdiseaseID_GDB = []
inh = []
mark = []
mim = []

with open(geneBaseFile,'r') as f:
    for line in f:
        line = line.replace('\n','')
        temp = line.split('\t')
        ConfidenceLevel.append(temp[0])
        TargetGroup.append(temp[1])
        geneBasedRef.append(temp[2])
        Disease.append(temp[3])
        SZAdiseaseID_GDB.append(temp[4])
        inh.append(temp[5])
        mark.append(temp[6])
        mim.append(temp[7])
geneBasedRefHeadings = [ConfidenceLevel[0],TargetGroup[0],geneBasedRef[0],Disease[0],inh[0],mark[0],mim[0]]

# Remove the first line for each col, which is the heading 
geneBasedRef.pop(0)
Disease.pop(0)
ConfidenceLevel.pop(0)
TargetGroup.pop(0)
SZAdiseaseID_GDB.pop(0)
inh.pop(0)
mim.pop(0)
mark.pop(0)

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

newContentHeadings = ['SZAID','Database','Database_type','SZAdiseaseID','ClinVar_CLNSIG','ClinVar_CLNDN','SZAreportCategory','ClinVar_ReviewStar','ClinVar_ID','selCSQ','Zygosity','Genotype']+geneBasedRefHeadings+headings;
newContent = '\t'.join(newContentHeadings)+'\n'


#%%Cell 4 report A is for Clinvar Variants
# Part A

for i in range(0,len(reportA)):
#for i in range(0,1):
    HGMD_flag = False
    if reportA[i].startswith('HGMD_'):
        HGMD_flag = True
        reportA[i] = reportA[i].replace('HGMD_','')

    #iText sepatated different transcript annotations
    iText = [s for s in reportA[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
    iREF = reportA[i].split('\t')[3].split(',')
    iALT = reportA[i].split('\t')[4].split(',')
    iGenoTypeList = iREF+iALT
    #check there is no missingness in the genotype field
    if "." not in reportA[i].split('\t')[9].split(':')[0]:
        iGenotypInd1 = int(reportA[i].split('\t')[9][0])
        iGenotypInd2 = int(reportA[i].split('\t')[9][2])
        #iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
        #8.21  for some male variants on chrX which only had one genotype
        if iGenotypInd2 <= len(iGenoTypeList)-1 :
            iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
        else:
            iGenotype = iREF[0]+'>'+iGenoTypeList[iGenotypInd1]
        iTemp = reportA[i].replace('\n','').split('\t')
        iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
        #sys.exit('!')
        if iGenotypInd1 == iGenotypInd2:
            iZygo = 'Homozygous'
        elif iGenotypInd1 == 0 or iGenotypInd2 == 0 or iGenotypInd2 > len(iGenoTypeList)-1:
            iZygo = 'Heterozygous'
        # Heterozygous should be consider if different variations occur on the same gene but different chromosomes (trans)
        # else:
        #     iZygo = 'Compound heterozygous'
        if check_male_flag == "Male":
            if "chrX" in reportA[i] and iZygo == 'Heterozygous':
                iZygo = "Hemizygous"
        for j in range(0,len(iText)):
            jText = iText[j].split('|')
            jGenes1 = jText[3]

            #jGenes2Temp = jText[41].split('&')
            jGenes2Temp = jText[colNames_CSQ.index("ClinVar_GENEINFO")].split('&')
            jGenes2 = []
            for k in range(0,len(jGenes2Temp)):
                kTemp = jGenes2Temp[k].split(':')
                jGenes2.append(kTemp[0])
            jGenes = list(set([jGenes1]+jGenes2))
            jGenes = [gene for gene in jGenes if gene != ""]
    ######### To clear GWAS and pharmaco ######
    ###The reason to do this: GWAS,pharmaco, Clinvar can be overlapped, so several SZAvar can be actually one, just showing in both clinvar and gwas or ..
    ###This will have issues like ValueError: 'ClinP_LP_var' is not in list
    ###Is because in one of the transcripts, the way of express alternative allele is different with what clinvar used.
    ## expression of "-" will annotated Clinvar_P_LP, while "TGCTGC"
            #print(jText[colNames_CSQ.index('Database_Type')].split('&'))
            
            #For CLN Review Star 
            Review = jText[colNames_CSQ.index('ClinVar_CLNREVSTAT')]
            if Review == 'criteria_provided&_conflicting_classifications':
                ReviewStar = 1
            elif Review == 'criteria_provided&_single_submitter':
                ReviewStar = 1
            elif Review == 'criteria_provided&_multiple_submitters&_no_conflicts':
                ReviewStar = 2
            elif Review == 'reviewed_by_expert_panel':
                ReviewStar = 3
            elif Review == 'practice_guideline':
                ReviewStar = 4
            else:
                ReviewStar = 0
            

            # For jCLNSIG
            if HGMD_flag == True:
                scoreFlag = 4 # same as the ClinVar Likely_pathogenic 
                jCLNSIG = "HGMD"
            
            elif HGMD_flag == False and jText[colNames_CSQ.index('Database_Type')]!= '':
                jCLINVARind = jText[colNames_CSQ.index('Database_Type')].split('&').index('ClinP_LP_var') # Find out the location index for ClinVar
                jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jCLINVARind]
                #jDatabase = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
                jCLNID = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
                jSZAvarID = jText[colNames_CSQ.index('Database_SZAID')].split('&')[jCLINVARind]
                #if jSZAvarID = 'SZAvar178179':
        ######### To clear GWAS and pharmaco #########
                jCLNSIG = jText[colNames_CSQ.index('ClinVar_CLNSIG')] #P/LP/Conflicting
                #print(jCLNSIG)
                #jCLNSIG = jText[colNames_CSQ.index('ClinVar_CLNSIG')].split('&')
                #print(jCLNSIG)
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ') # ClinVar Disease name
                if 'Pathogenic' in jCLNSIG:
                    scoreFlag = 5
                elif 'Likely_pathogenic' in jCLNSIG:
                    scoreFlag = 4
                #elif 'risk_factor' in jCLNSIG:
                #    scoreFlag = 1
                elif 'Conflicting_interpretations_of_pathogenicity' in jCLNSIG:
                    scoreFlag = 1
                elif 'Conflicting_classifications_of_pathogenicity' in jCLNSIG:
                    pathogenicity = jText[colNames_CSQ.index('ClinVar_CLNSIGCONF')].split('&')
                    pathogenicity = [re.sub(r'\(\d+\)', '', s) for s in pathogenicity]
                    if 'Pathogenic' in pathogenicity or 'Likely_pathogenic' in pathogenicity:    
                        scoreFlag = 2
                    else:
                        scoreFlag = 1
                else:
                    scoreFlag = 0
                #if 'KRT75' in jGenes:
                #    sys.exit('!')
                # Extract the ClinVar genes that also exist in GeneDB in each line 
            else:
                sys.exit('!')

            ## 1. gene is in GeneDB
            if scoreFlag > 0 and any([x for x in jGenes if x in geneBasedRef]): 
                jGeneList = [x for x in jGenes if x in geneBasedRef] # store all associated genes as a list 
                jInds = []
                #if 'HBG2' in jGeneList:
                #    sys.exit('!')
                for gene in range(0,len(jGeneList)):
                    kTemp = jGeneList[gene]
                    jInds = jInds +[ii for ii, x in enumerate(geneBasedRef) if x == kTemp] # Find gene location in GeneDB
                    # if gene in geneBasedRef(GeneDB) is same as kTemp (in jGeneList, which is the ClinVar gene),
                    # then return the index (location) as ii, and store all locations in jInds

                ######### test

                # Extract the ClinVar disease in each line 
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ')
                kClinVarDiseases = jClinVar_CLNDN.split('|') # All the ClinVar disease in each line 
                #print(kClinVarDiseases)
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                for c in range(0,len(kClinVarDiseases)):
                    kClinVarDisease_temp = kClinVarDiseases[c].replace(',','')
                    if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1: 
                        kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp == x.replace(',','')]
                        #kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp in x.replace(',','') or x.replace(',','') in kClinVarDisease_temp]
                        for element in kIndex:
                            jVarSZAdiseaseIDs.append(DiseaseID_DSDB[element])
                            jVarDiseaseNames.append(DiseaseName_DSDB[element])

                   #sys.exit('!')
                ## problem here!!! (should match Clinvar disease type?)
                ## fixed by if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs: SZAdiseaseID_GDB is the disease from GeneDB, to see
                # whether it's matching with the ClinVar_CLNDN disease jVarSZAdiseaseIDs
                # Check if the ClinVar gene-disease pair is also in the GeneDB, and assign the confidence levels 
                for k in range(0,len(jInds)):
                    if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs: # if the GeneDB disease is in DiseaseDB  disease (exactly match)
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
                        
                        if HGMD_flag == False:
                            jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                            jCLNID = jText[colNames_CSQ.index('Database')]
                            newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                        Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n'
                        elif HGMD_flag == True:
                            newContent = newContent+'\t'.join(['NA','NA','HGMD',jSZADiseaseID,jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n'  
                        
                        # Diseaes in newContent comes from the GeneDB
                                        
                    else:
                        GeneDB_szaID = SZAdiseaseID_GDB[jInds[k]]
                        GeneDB_Disease = Disease[jInds[k]]
                        GeneDB_Disease_word = set(GeneDB_Disease.upper().split())

                        if len(jVarDiseaseNames) > 0:
                            for name in range(len(jVarDiseaseNames)):
                                clinvar_disease_word = set(jVarDiseaseNames[name].upper().split())
                                matches = False

                                if clinvar_disease_word.issubset(GeneDB_Disease_word):
                                    matches = True
                                elif GeneDB_Disease_word.issubset(clinvar_disease_word):
                                    matches = True
                            
                                if matches == True:
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

                                    if HGMD_flag == False:
                                        jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                                        jCLNID = jText[colNames_CSQ.index('Database')]                                    
                                        newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                                    Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n'
                                    elif HGMD_flag == True:
                                        newContent = newContent+'\t'.join(['NA','NA','HGMD',jSZADiseaseID,jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                    Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n' 

                                else:
                                ########here, if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs is false, then no SZAscore value
                                    SZAscore = scoreFlag
                                    jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]

                                    if HGMD_flag == False:
                                        jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                                        jCLNID = jText[colNames_CSQ.index('Database')]
                                        newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar),jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                                        Disease[jInds[k]],inh[jInds[k]],mark[jInds[k]], mim[jInds[k]]]+ iTemp)+'\n'
                                    elif HGMD_flag == True:
                                        newContent = newContent+'\t'.join(['NA','NA','HGMD',jSZADiseaseID,jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                                    Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n'                             
                                    
                                #Diseaes in newContent comes from the GeneDB
                        else: # for variants do not have ClinVar disease
                            SZAscore = scoreFlag 
                            if HGMD_flag == True:
                                newContent = newContent+'\t'.join(['NA','NA','HGMD',SZAdiseaseID_GDB[jInds[k]],jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n'                           

                            else:
                                newContent = newContent+'\t'.join(['NA',jCLNID,jDatabaseType,SZAdiseaseID_GDB[jInds[k]],jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID, iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],
                                Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]] +iTemp)+'\n' 


                            
            elif scoreFlag >0 and any(jGenes): # If there is jGenes, but not included in GeneDB
                #sys.exit('Found a gene set with '+', '.join(jGenes)+' that is not included in GeneDB ')
                SZAscore = scoreFlag
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ')
                kClinVarDiseases = jClinVar_CLNDN.split('|')
                #print(kClinVarDiseases)
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                jSZADiseaseID = []
                jDisease = []
                for k in range(0,len(kClinVarDiseases)):
                    #kClinVarDisease_temp = kClinVarDiseases[k]
                    kClinVarDisease_temp = kClinVarDiseases[k].replace(",","")
                    if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1:
                        kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]
                        jSZADiseaseID_tmp = DiseaseID_DSDB[kIndex[0]]
                        jSZADiseaseID.append(jSZADiseaseID_tmp)
                        jDisease_tmp = DiseaseName_DSDB[kIndex[0]]
                        jDisease.append(jDisease_tmp)
                    else:
                        jSZADiseaseID_tmp = ''
                        jSZADiseaseID.extend(jSZADiseaseID_tmp)
                        jDisease_tmp = ''
                        jDisease.append(jDisease_tmp)
                jDisease = [x for x in jDisease if x != '']
                jSZADiseaseID = [x for x in jSZADiseaseID if x != '']
                if len(jDisease) == 1: 
                    jDisease = ''.join(jDisease)
                    jSZADiseaseID = ''.join(jSZADiseaseID)

                    if HGMD_flag == False:
                        jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                        jCLNID = jText[colNames_CSQ.index('Database')]
                        newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease,"None","None","None"]+iTemp)+'\n'
                    elif HGMD_flag == True:
                        newContent = newContent+'\t'.join(['NA','NA','HGMD',jSZADiseaseID,jCLNSIG,'NA',str(SZAscore),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease,"None","None","None"]+iTemp)+'\n'
                
                
                elif len(jDisease) > 1:
                    for m in range(len(jDisease)):
                        mjDisease = ''.join(jDisease[m])
                        mjSZADiseaseID = ''.join(jSZADiseaseID[m])
                        if HGMD_flag == False:
                            jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                            jCLNID = jText[colNames_CSQ.index('Database')]
                            newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,mjSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),mjDisease,"None","None","None"]+iTemp)+'\n'
                        elif HGMD_flag == True:
                            newContent = newContent+'\t'.join(['NA','NA','HGMD',mjSZADiseaseID,jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),mjDisease,"None","None","None"]+iTemp)+'\n'
                else: # len(jDisease) == 0, for HGMD variant without ClinVar disease 
                    if HGMD_flag == False:
                        newContent = newContent+'\t'.join(['NA',jCLNID,jDatabaseType,"NA",jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),"None","None","None","None"]+iTemp)+'\n'

                    else:
                        newContent = newContent+'\t'.join(['NA','NA','HGMD',"NA",jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),"None","None","None","None"]+iTemp)+'\n'
                    

            elif scoreFlag >0:
                SZAscore = scoreFlag
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ')
                kClinVarDiseases = jClinVar_CLNDN.split('|')
                jVarSZAdiseaseIDs = []
                jVarDiseaseNames = []
                jSZADiseaseID = []
                jDisease = []
                for k in range(0,len(kClinVarDiseases)):
                    kClinVarDisease_temp = kClinVarDiseases[k].replace(",","")
                    if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1:
                        kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]
                        #jSZADiseaseID = DiseaseID_DSDB[kIndex[0]]
                        #jDisease = DiseaseName_DSDB[kIndex[0]]
                        jSZADiseaseID_tmp = DiseaseID_DSDB[kIndex[0]]
                        jSZADiseaseID.append(jSZADiseaseID_tmp)
                        jDisease_tmp = DiseaseName_DSDB[kIndex[0]]
                        jDisease.append(jDisease_tmp)
                    else:
                        jSZADiseaseID_tmp = ''
                        #jDisease = jClinVar_CLNDN
                        jSZADiseaseID.extend(jSZADiseaseID_tmp)
                        jDisease_tmp = ''
                        jDisease.append(jDisease_tmp)
                        #sys.exit('!')
                jDisease = [x for x in jDisease if x != '']
                jSZADiseaseID = [x for x in jSZADiseaseID if x != '']
                #sys.exit('!')
                if len(jDisease) == 1: 
                    jDisease = ''.join(jDisease)
                    jSZADiseaseID = ''.join(jSZADiseaseID)
                    if HGMD_flag == False:
                        jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                        jCLNID = jText[colNames_CSQ.index('Database')]
                        newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease,"None","None","None"]+iTemp)+'\n'
                    elif HGMD_flag == True:
                        newContent = newContent+'\t'.join(['NA','NA','HGMD',jSZADiseaseID,jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease,"None","None","None"]+iTemp)+'\n'
                
                elif len(jDisease) > 1:
                    for m in range(len(jDisease)):
                        mjDisease = ''.join(jDisease[m])
                        mjSZADiseaseID = ''.join(jSZADiseaseID[m])
                        if HGMD_flag == False:
                            jDatabaseType = jText[colNames_CSQ.index('Database_Type')]
                            jCLNID = jText[colNames_CSQ.index('Database')]
                            newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,mjSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),mjDisease,"None","None","None"]+iTemp)+'\n'
                        elif HGMD_flag == True:
                            newContent = newContent+'\t'.join(['NA','NA','HGMD',mjSZADiseaseID,jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),mjDisease,"None","None","None"]+iTemp)+'\n'
                else:
                    if HGMD_flag == True:
                        newContent = newContent+'\t'.join(['NA','NA','HGMD',"NA",jCLNSIG,'NA',str(0),str(ReviewStar), 'NA',iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),"None","None","None","None"]+iTemp)+'\n'
                    else:
                        newContent = newContent+'\t'.join(['NA',jCLNID,jDatabaseType,"NA",jCLNSIG,jClinVar_CLNDN,str(SZAscore),str(ReviewStar), jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),"None","None","None","None"]+iTemp)+'\n'


## Compound Heterozygous analyses (ClinVar and HGMD variants will both count)
clinvar_headings = newContent.split('\n')[0].split('\t')
line = [clinvar.split('\t') for clinvar in newContent.split('\n')[1:]]
clinvar_df = pd.DataFrame(data=line, columns=clinvar_headings)
clinvar_df["Compound_heterozygous"] = "No"
compound_tmp_df = clinvar_df[['Genes','#CHROM','POS','REF','ALT',clinvar_df.columns[28]]].drop_duplicates().dropna()
compound_tmp_df.iloc[:, 5] = compound_tmp_df.iloc[:, 5].apply(lambda x: x.split(':')[0])
compound_tmp_df['Genotype'] = None
compound_tmp_df = compound_tmp_df.reset_index(drop=True)

#Extract the variation and zygosity for compound heterozygous 
for i in range(len(compound_tmp_df)):
    genotype0 = compound_tmp_df.iloc[i,5][0]
    if len(compound_tmp_df.iloc[i,5]) > 1:
        genotype1 = compound_tmp_df.iloc[i,5][1]
        genotype2 = compound_tmp_df.iloc[i,5][2]
    else:
        genotype1 = genotype0
        genotype2 = genotype0
    if genotype0 != genotype2:
        if genotype1 == "/":
            compound_tmp_df.at[i,"Genotype"] = "Heterozygous(Unphased)"
        elif genotype1 == "|":
            compound_tmp_df.at[i,"Genotype"]  = "Heterozygous(Phased)"

compound_tmp_df = compound_tmp_df.dropna()
grouped = compound_tmp_df.groupby('Genes').filter(lambda x: len(x) > 1)

unique_genes = grouped['Genes'].unique()
for gene in unique_genes:
    tmp_df = grouped[grouped['Genes'] == gene]
    genotype = tmp_df[tmp_df['Genes'] == gene]['Genotype']
    if 'Heterozygous(Unphased)' in genotype.values:
        for index, row in tmp_df.iterrows():
            ref_tmp = row["REF"]
            chrom_tmp = row['#CHROM']
            alt_tmp = row['ALT']
            pos_tmp = row['POS']
            index = clinvar_df[(clinvar_df['Genes'] == gene) & 
                            (clinvar_df['REF'] == ref_tmp) &
                            (clinvar_df['#CHROM'] == chrom_tmp) &
                            (clinvar_df['POS'] == pos_tmp) &
                            (clinvar_df['ALT'] == alt_tmp)].index
            if len(index) > 0:
                clinvar_df.loc[index, 'Compound_heterozygous'] = 'Uncertained ' + row.iloc[5]
    else:
        clinvar_df.loc[index, 'Compound_heterozygous'] = 'Yes'


## Change df into NewContent
newContent = '\t'.join(newContentHeadings + ["Compound_heterozygous"])+'\n'
for index, row in clinvar_df.iterrows():
    newContent += '\t'.join(row.astype(str)) + '\n'
                        

#%%Cell 5 report C GWAS

if run_GWAS== True:
    from datetime import datetime
    newLineHeadings = ['Trait','Genotype','Zygosity','risk_genotype','SZAvarID','GWASID','OR','Pval','PUBMED','Study accession ID','MAPPED_GENE',"CONTEXT"];
    newLine = '\t'.join(newLineHeadings)+'\n'
    with open(outFile.replace('.txt','_GWAS.txt'),'w') as f:
        f.write(newLine)
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportC)) + ' lines of GWAS variants' )
    # GWAS part
    for i in range(0,len(reportC)):
        #print(i)
        current_time = now.strftime("%H:%M:%S")
        iText = [s for s in reportC[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
        iREF = reportC[i].split('\t')[3].split(',')
        iALT = reportC[i].split('\t')[4].split(',')
        iGenoTypeList = iREF+iALT
        #to fix the problem in swedish vcf files:
        #if "." not in reportC[i].split('\t')[9]:

        ####For turkish cohort:
        if "." not in reportC[i].split('\t')[9].split(':')[0]:
            iGenotypInd1 = int(reportC[i].split('\t')[9][0])
            iGenotypInd2 = int(reportC[i].split('\t')[9][2])
            #some variants have genotype 1/2 but with ALT only one, strange...
            #iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            try:
                iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            except:
                iGenotype = iREF[0]+'>'+iALT[0]+':'+reportC[i].split('\t')[9][0]+reportC[i].split('\t')[9][1]+reportC[i].split('\t')[9][2]
            iTemp = reportC[i].replace('\n','').split('\t')
            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
            # Get the GWAS ID and SZAvarID from the first sub transcript, named jText
            ## here we seems not care about genes, transcripts, only focuse on variant so I deleted for loop
            jText = iText[0].split('|')
            if jText[colNames_CSQ.index('Database_Type')]!= '':
                jGWASind = jText[colNames_CSQ.index('Database_Type')].split('&').index('GWAS_var')
                jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jGWASind]
                #jDatabase = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
                jGWASID = jText[colNames_CSQ.index('Database')].split('&')[jGWASind]
                jSZAvarID = jText[colNames_CSQ.index('Database_SZAID')].split('&')[jGWASind]
                ######### To match the GWAS_DB
                #one sza variant can mapped to different traits and records
                ##first: explor how many records matched for one szaid but with different study(so they have different p value, OR):
                # Filter GWAS data for the SZAvarID
                GWAS_db_filtered = GWAS_db[GWAS_db['SZAID'] == jSZAvarID]
                #jStudy_accession_list = GWAS_db[GWAS_db['SZAID']== jSZAvarID]['STUDY.ACCESSION'].to_list()
                for k, row in GWAS_db_filtered.iterrows():
                    jOR = row['OR.or.BETA']
                    jPval = row['P.VALUE']
                    jrisk_genotype = row['risk_genotype']
                    jPUBMED = row['PUBMEDID']
                    jTrait = row['DISEASE.TRAIT']
                    jMappedGene = row['MAPPED_GENE']
                    jVarType = row['CONTEXT']
                    jStudy_accession = row['STUDY.ACCESSION']
                    #only keep p value <= 5e-8 for report
                    if iGenotype == jrisk_genotype and jPval <= 5e-8:
                        newLine_temp = '\t'.join([jTrait, iGenotype, iZygo, jrisk_genotype, jSZAvarID, jGWASID, str(jOR), str(jPval), str(jPUBMED), str(jStudy_accession), jMappedGene, jVarType]) + '\n'
                        with open(outFile.replace('.txt', '_GWAS.txt'), 'a') as f:
                            f.write(newLine_temp)



        if i%1000 == 0:
            print(str(i)+' lines processed!')

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



#%% Cell 6 report D for pharmaco-annotation

#if run_Pharmaco == True:

#    from datetime import datetime
#    now = datetime.now()
#    newLineHeadings = ['Drug.s.','Genotype','Zygosity','pharma_genotype','SZAvarID','Annotation.Text','Score','Phenotype.Category','Phenotype.s.','Gene','Level.of.Evidence','Level.Modifiers','URL'];
#    newLine = '\t'.join(newLineHeadings)+'\n'
#    with open(outFile.replace('.txt','_pharmaco.txt'),'w') as f:
#        f.write(newLine)
#    current_time = now.strftime("%H:%M:%S")
#    print("Current Time =", current_time)
#    print('processing ' + str(len(reportD)) + ' lines of pharmaco variants' )
#    # Pharma part
#    for i in range(0,len(reportD)):
#        current_time = now.strftime("%H:%M:%S")
#       iText = [s for s in reportD[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
#        iREF = reportD[i].split('\t')[3].split(',')
#        iALT = reportD[i].split('\t')[4].split(',')
#        iGenoTypeList = iREF+iALT
        ##some variants in swedish cohort has: ./0 like this chr7    87531302        rs2032582
        ## but on the previous line exactly same variants ,so make just skip the line
        #So that int() will gives error
        ### To fix this:
        ###only the normal variants with normal genotypes goes into outfile
        ##For Turkish
#       if "." not in reportD[i].split('\t')[9].split(':')[0]:
#            iGenotypInd1 = int(reportD[i].split('\t')[9][0])
#            iGenotypInd2 = int(reportD[i].split('\t')[9][2])
            #iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
#            try:
#                iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
#            except:
#               iGenotype = iREF[0]+'>'+iALT[0]+':'+reportD[i].split('\t')[9][0]+reportD[i].split('\t')[9][1]+reportD[i].split('\t')[9][2]
#            iTemp = reportD[i].replace('\n','').split('\t')
#            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
#            #sys.exit('!')
#            if iGenotypInd1 == iGenotypInd2:
#                iZygo = 'Homozygous'
#            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
#                iZygo = 'Heterozygous'
#            else:
#                iZygo = 'Compound heterozygous'
#                #print(iZygo)
            ##get the Pharma and sza id from sub transcripts, named jText
            ## here we seems not care about genes, transcripts, only focuse on variant so I deleted for loop
#            jText = iText[0].split('|')
#            if jText[colNames_CSQ.index('Database_Type')]!= '':
#                jPharmaind = jText[colNames_CSQ.index('Database_Type')].split('&').index('Pharma_var')
#                jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jPharmaind]
                #jDatabase = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
 #               jPharmaID = jText[colNames_CSQ.index('Database')].split('&')[jPharmaind]
 #               jSZAvarID = jText[colNames_CSQ.index('Database_SZAID')].split('&')[jPharmaind]
                ######### To match the Pharma_DB
                ### The pharma_genotype for one SZAID is multiple!!
                # so need to use genotype_ind
#              jpharma_genotype_list = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['pharma_genotype'].to_list()
#                if iGenotype in jpharma_genotype_list:
#                    genotype_ind = jpharma_genotype_list.index(iGenotype)
#                    jpharma_genotype = jpharma_genotype_list[genotype_ind]
#                    jDrug= Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Drug.s.'].to_list()[genotype_ind]
#                    jAnnotation = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Annotation.Text'].to_list()[genotype_ind]
#                    jScore = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Score'].to_list()[genotype_ind]
#                    jPheno_category = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Phenotype.Category'].to_list()[genotype_ind]
#                    jPheno = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Phenotype.s.'].to_list()[genotype_ind]
#                    jGene = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Gene'].to_list()[genotype_ind]
#                    jevidence = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Level.of.Evidence'].to_list()[genotype_ind]
#                    jlevel_modifier = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Level.Modifiers'].to_list()[genotype_ind]
#                    jURL = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['URL'].to_list()[genotype_ind]
#                    newLine_temp ='\t'.join([jDrug,iGenotype,iZygo,jpharma_genotype,jSZAvarID,jAnnotation,str(jScore),jPheno_category,jPheno,jGene,jevidence,jlevel_modifier,jURL])+'\n'
#                    with open(outFile.replace('.txt','_pharmaco.txt'),'a') as f:
#                        f.write(newLine_temp)
        #newLineHeadings = ['Drug(s)','Genotype','pharma_genotype','SZAvarID','Annotation Text','Score','Phenotype Category','Phenotype(s)','Gene','Level of Evidence','Level Modifiers','URL'];

#    now = datetime.now()
#    current_time = now.strftime("%H:%M:%S")
#    print("Current Time =", current_time)

#%% Cell 6 (new) report D for pharmaco-annotation

if run_Pharmaco == True:

    from datetime import datetime
    now = datetime.now()
    #newLineHeadings = ['Drug.s.','Genotype','Zygosity','pharma_genotype','SZAvarID','Annotation.Text','Score','Phenotype.Category','Phenotype.s.','Gene','Level.of.Evidence','Level.Modifiers','URL'];
    #newLine = '\t'.join(newLineHeadings)+'\n'
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportD)) + ' lines of pharmaco variants' )
    #Pharma_db['Genotype.Allele'] = Pharma_db['Genotype.Allele'].astype(str)
    #Pharma_db['Zygous'] = Pharma_db['Zygous'].astype(str)
    # Pharma part
    for i in range(0,len(reportD)):
    #for i in range(726,727):
        current_time = now.strftime("%H:%M:%S")
        iText = [s for s in reportD[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
        iREF = reportD[i].split('\t')[3].split(',')
        iALT = reportD[i].split('\t')[4].split(',')
        iGenoTypeList = iREF+iALT
        ##some variants in swedish cohort has: ./0 like this chr7    87531302        rs2032582
        ## but on the previous line exactly same variants ,so make just skip the line
        #So that int() will gives error
        ### To fix this:
        ###only the normal variants with normal genotypes goes into outfile
        ##For Turkish
        if "." not in reportD[i].split('\t')[9].split(':')[0]:
            iGenotypInd1 = int(reportD[i].split('\t')[9][0])
            iGenotypInd2 = int(reportD[i].split('\t')[9][2])
            genotype1 = iGenoTypeList[iGenotypInd1] + iGenoTypeList[iGenotypInd2]
            genotype_sorted = ''.join(sorted(genotype1))
            #iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            try:
                iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            except:
                iGenotype = iREF[0]+'>'+iALT[0]+':'+reportD[i].split('\t')[9][0]+reportD[i].split('\t')[9][1]+reportD[i].split('\t')[9][2]
            iTemp = reportD[i].replace('\n','').split('\t')
            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
            #sys.exit('!')
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
                #print(iZygo)
            ##get the Pharma and sza id from sub transcripts, named jText
            ## here we seems not care about genes, transcripts, only focuse on variant so I deleted for loop
            jText = iText[0].split('|')
            if jText[colNames_CSQ.index('Database_Type')]!= '':
                content = jText[colNames_CSQ.index('Database_Type')].split('&')
                if 'Pharma_var' in content:
                    jPharmaind = jText[colNames_CSQ.index('Database_Type')].split('&').index('Pharma_var')
                    jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jPharmaind]
                    #jDatabase = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
                    jPharmaID = jText[colNames_CSQ.index('Database')].split('&')[jPharmaind]
                else:
                    IDs = jText[colNames_CSQ.index('Existing_variation')].split('&')
                    index = next((i for i, ID in enumerate(IDs) if ID.startswith('rs')), None)
                    if index != None:
                        jPharmaID = IDs[index]
                    else:
                        jPharmaID = 'na'
            else:
                IDs = jText[colNames_CSQ.index('Existing_variation')].split('&')
                if IDs != ['']:
                    index = next((i for i, ID in enumerate(IDs) if ID.startswith('rs')), None)
                    if index != None:
                        jPharmaID = IDs[index]
                    else: 
                        jPharmaID = 'na'
                else:
                    jPharmaID = 'na'
                
            jPharmaID = str(jPharmaID)


            ichr = iTemp[0]
            ipos = iTemp[1]
            iref = iTemp[3]
            ialt = iTemp[4]
            ivariation = f"{ichr}_{ipos}_{iref}_{ialt}"


            for r in range(len(Pharma_db)):
                tmp_id = Pharma_db.loc[Pharma_db.index[r], 'rsID']
                if jPharmaID == tmp_id:
                    if '*' in Pharma_db.loc[Pharma_db.index[r], 'Genotype.Allele']:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = iGenoTypeList[iGenotypInd1] + '/' + iGenoTypeList[iGenotypInd2]
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo
                    else:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = genotype_sorted
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo
                if ivariation == tmp_id:
                    if '*' in Pharma_db.loc[Pharma_db.index[r], 'Genotype.Allele']:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = iGenoTypeList[iGenotypInd1] + '/' + iGenoTypeList[iGenotypInd2]
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo
                    else:
                        Pharma_db.loc[Pharma_db.index[r],"Genotype"] = genotype_sorted
                        Pharma_db.loc[Pharma_db.index[r],"Zygous"]= iZygo

    def filter_genotypes(group):
        return not (group['Genotype'] == 'a').any()
    Pharma_db_final = Pharma_db.groupby('Variant.Haplotypes').filter(filter_genotypes)

    Pharma_db_final = Pharma_db_final[
    ~((Pharma_db_final['Variant.Haplotypes'].str.startswith('rs')) & 
      (Pharma_db_final['Genotype'] != Pharma_db_final['Genotype.Allele']))]
                
    pgx_filename =outFile.replace('.txt','_PGx.txt') 

    Pharma_db_final.to_csv(pgx_filename, index=False,sep='\t')
  

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

#%% Cell 7 Traits report
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
    irs_id_set = {s.split('|')[colNames_CSQ.index('Existing_variation')].split('&')[0] for s in reportE}
    for i in range(0,len(reportE)):
        current_time = now.strftime("%H:%M:%S")
        ##iText means all transcripts
        iText = [s for s in reportE[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
        irs_id = iText[0].split("|")[colNames_CSQ.index('Existing_variation')].split('&')[0]
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
                #print(iZygo)
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
            print(row['variants'] + "is a wild type genotype")
            newLine_temp = '\t'.join([row['Traits name'], row['category'], row['genes'], row['variants'], row['Description'], row['REF']+'/'+row['REF'], row['Genotype_Description'], 'homo_ref']) + '\n'
            with open(outFile.replace('.txt', '_traits.txt'), 'a') as f:
                f.write(newLine_temp)
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
        local_file_path = outFile.replace('.txt', '_traits.txt')
        remote_file_path = '/home/ubuntu/bckndsrv/RawDataOutputs(GenomicSZA)/' + local_file_path.split('/')[-1]

        # Transfer the file
        scp.put(local_file_path, remote_file_path)
        print("File transferred successfully to AWS directory.")

        # Close SCP and SSH sessions
        scp.close()
        ssh.close()
#%% Cell 8 Part B processing and generate outFile
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

varCount = 0
for i in range(0,len(reportB)):
    iText = [s for s in reportB[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
    iREF = reportB[i].split('\t')[3].split(',')
    iALT = reportB[i].split('\t')[4].split(',')
    iGenoTypeList = iREF+iALT
    ####For swedish cohort
    #if "." not in reportB[i].split('\t')[9]:
    #but it's not gonna work for turkish format
    ####For turkish cohort:
    if "." not in reportB[i].split('\t')[9].split(':')[0]:
        iGenotypInd1 = int(reportB[i].split('\t')[9][0])
        iGenotypInd2 = int(reportB[i].split('\t')[9][2])
        try:
            iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
        except:
            iGenotype = iREF[0]+'>'+iALT[0]+':'+reportB[i].split('\t')[9][0]+reportB[i].split('\t')[9][1]+reportB[i].split('\t')[9][2]
        if iGenotypInd1 == iGenotypInd2:
            iZygo = 'Homozygous'
        elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
            iZygo = 'Heterozygous'
        else:
            iZygo = 'Compound heterozygous'
        ##New content replace the Zygosity of the chrX for males
        if check_male_flag == "Male":
            if "chrX" in reportB[i] and iZygo == 'Heterozygous':
                iZygo = "Hemizygous"
        varFlag = 0
        iTemp = reportB[i].replace('\n','').split('\t')
        iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
        for j in range(0,len(iText)):
            jText = iText[j].split('|')
            if (jText[colNames_CSQ.index('MAX_AF')] == '' or float(jText[colNames_CSQ.index('MAX_AF')])<0.05):
                ##Update: to clear the logics here need to add brackets
                ## double check
                ## (20å24-5-29) need to change the context into:
                    # 1. for Loss of Functions: High IMPACT and HC_LoF(LOFTEE) (LOFTEE can be filtered later in R)
                    # 2. for Missense variants: AM_likely_pathogenic and REVEL > 0.75
                    # 3. for Splicing variants: SpliceAI > 0.8 
                    # 4. Overall scores: Bayels_addAF_score > 0.0692655 and BayesDel_noAF_score > -0.0570105
                
                if (jText[colNames_CSQ.index('IMPACT')] == 'HIGH' or (
                        jText[colNames_CSQ.index('ada_score')] != '' and float(jText[colNames_CSQ.index('ada_score')])>0.6) or (
                        jText[colNames_CSQ.index('rf_score')] != '' and float(jText[colNames_CSQ.index('rf_score')])>0.6) or (
                        jText[colNames_CSQ.index('REVEL')]!= '' and float(jText[colNames_CSQ.index('REVEL')])>0.75) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')])>0.5) or (
                        jText[colNames_CSQ.index('BayesDel_addAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_addAF_score')])>0.0692655) or (
                        jText[colNames_CSQ.index('BayesDel_noAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_noAF_score')])>-0.0570105) or (
                        jText[colNames_CSQ.index('am_class')]== 'likely_pathogenic' and float(jText[colNames_CSQ.index('am_pathogenicity')])> 0.564) or (
                        jText[colNames_CSQ.index('LoF')]== 'HC')):

                    jGenes = jText[3].split('&')
                    if any([x for x in jGenes if x in geneBasedRef]):
                        varFlag = 1
                        jGeneList = [x for x in jGenes if x in geneBasedRef]
                        #if 'KLC1' == jGeneList[0]:
                        #    sys.exit('!')
                        jInds = []
                        for k in range(0,len(jGeneList)):
                            kTemp = jGeneList[k]
                            jInds = jInds +[i for i, x in enumerate(geneBasedRef) if x == kTemp]
                        for k in range(0,len(jInds)):
                            scoreFlag = 2+(jText[colNames_CSQ.index('IMPACT')]=='HIGH')*1
                            #SZAscore =0
                            if ConfidenceLevel[jInds[k]] == 'High_confidence':
                                SZAscore = 10+scoreFlag
                            elif ConfidenceLevel[jInds[k]] == 'Moderate_confidence':
                                SZAscore = 5+scoreFlag
                            elif ConfidenceLevel[jInds[k]] == 'Low_confidence':
                                SZAscore = scoreFlag
                            else:
                                print(jGenes)
                                sys.exit('Unknown confidence score level for the report!')
                            jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                            #print(jGenes)
                            #SZAscore = 2+('High_confidence' == ConfidenceLevel[jInds[k]])*5+(jText[colNames_CSQ.index('IMPACT')]=='HIGH')*1
                            newContent = newContent+'\t'.join(
                                ['NovelTrans'+headings[9]+'_'+str(varCount+1),'NA','NA',
                                jSZADiseaseID,'NA','NA',str(SZAscore),'NA','NA',iText[j],iZygo,
                                iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],
                                geneBasedRef[jInds[k]],Disease[jInds[k]],inh[jInds[k]], mark[jInds[k]], mim[jInds[k]]]+iTemp + ["No"])+'\n'
                                # The last "No" is for compound_heterozygous
                
                
    if varFlag == 1:
        varCount = varCount+1

print("putative variants matched the GeneDB:",varCount)

# write final report
outfid = open(outFile,'w')
outfid.write(newContent)
outfid.close()

print("outFile generated. ")


#%%(when do batch processing)
#Get the original sample outfile names from the python_output directory
#file_count = 0
#names =[]
#for name in os.listdir("/Users/xiyas/Report_genome/original_outfile"):
#for name in os.listdir("/Users/xiyas/Report_genome/Sample_data_SZA"):
    #print(name)
    #file_count +=1
    #print(file_count)
#    if "outfile.txt" in name:
#        name = os.path.join("/Users/xiyas/Report_genome/original_outfile", name)
        #name = os.path.join("/proj/snic2020-16-69/nobackup/WGS_SZA/python_outputfile/outfile_dir", name)
#        names.append(name)

#%%Cell 10 separate the selCSQ generate a final report "final_report_sp"with all csq informations
#for k in names:
# outFile = outFile
# print("processing"+outFile)

# ##################
outFile_sp = outFile.replace('.txt','_sp.txt')
# read original final report
final_report_sp  = pd.read_csv(outFile, sep="\t")
# split the selCSQ
#final_report_sp[colNames_CSQ]=final_report_sp['selCSQ'].apply(lambda x: pd.Series(str(x).split("|")))
# save as new format
split_df = final_report_sp['selCSQ'].str.split("|", expand=True)
split_df.columns = colNames_CSQ
final_report_sp = pd.concat([final_report_sp, split_df], axis=1)
final_report_sp.to_csv(outFile_sp, sep = "\t",index = None, header=True,)
#Cell 5 Adding Inheritance and disease info
## adding inheritance from OMIM
outFile_sp_Inheritance = outFile.replace('.txt','_sp_Inheritance.txt')
NewHeadings = []
NewHeadings = list(final_report_sp.columns) + ['Inheritances'] + ['DiseaseInfo']

with open((outFile_sp_Inheritance),'w') as f2:
    f2.write('\t'.join(NewHeadings)+'\n')

#matching the inheritance
with open(outFile_sp,'r') as f1:
    next(f1)
    #linecount= len(f.readlines())
    for line in f1:
        line = line.replace('\n','')
        temp = line.split('\t')
        iSZAdiseaseID = temp[3]
        #print(temp[3])
        if 'SZA' not in temp[3]:
            continue
        #iGene = temp[13]
        iGene = temp[31] # Change the variable from 'Genes' to 'SYMBOL' 
        iOMIM_num = DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseMIM'].to_list()[0]
        #print(iOMIM_num)
        if iOMIM_num:
            iheritance = OMIM_Inheritance_DB[(OMIM_Inheritance_DB['phenotypeMimNumber']== iOMIM_num) &
                                        (OMIM_Inheritance_DB['approvedGeneSymbol']== iGene)]['inheritances'].to_list()
            if iheritance:
                iheritance = iheritance[0]
            else:
                iheritance = "No matched diseaseMIM in OMIM/OMIM not provided Inheritance"
        #iDiseaseInfo writing" now only use the column named OMIM_Description in the disease database
        #iDiseaseInfo = "Disease description(MONDO):" +DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['Mondo_Description'].to_list()[0]+ '|' + "Disease description(OMIM):" +DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['OMIM_Description'].to_list()[0]+ '|' + "Reference:"+ "DiseaseName:" + DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseName'].to_list()[0]+ '|' + "DiseaseSource:" + DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['SourceName'].to_list()[0]+'|' +"Disease SourceID:" +DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['SourceID'].to_list()[0]+'|' + "OMIM number:"+ DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseMIM'].to_list()[0]
        iDiseaseInfo = DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['OMIM_Description'].to_list()[0]+ '|' + "Reference:"+ "DiseaseName:" + DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseName'].to_list()[0]+ '|' + "DiseaseSource:" + DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['SourceName'].to_list()[0]+'|' +"Disease SourceID:" +DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['SourceID'].to_list()[0]+'|' + "OMIM number:"+ DiseaseDB[DiseaseDB['SZAdiseaseID']== iSZAdiseaseID]['DiseaseMIM'].to_list()[0]
        line = line + "\t" + iheritance +"\t"+iDiseaseInfo
        with open((outFile_sp_Inheritance),'a') as f2:
            f2.write(line+ '\n')

#Cell 6 split the 3 reports based on scores before removing duplicated records
### 8.7 add a nodup4 file that contains all the matching records
final_report_sp_Inheritance  = pd.read_csv(outFile_sp_Inheritance, sep="\t",dtype=str)
#final_report_general_1 = final_report_sp_Inheritance[final_report_sp_Inheritance['SZAreportCategory'].isin(['14','15'])]
#    final_report_general_1 = final_report_general_1[["Target.group","Disease","ClinVar_CLNDN","Genes","SZAID","Existing_variation","VARIANT_CLASS","Consequence","ClinVar_CLNHGVS","#CHROM","POS","ClinVar","REF","ALT","ClinVar_CLNSIG","Inheritances","DiseaseInfo","Zygosity","ClinVar_GENEINFO"]]
#final_report_general_2 = final_report_sp_Inheritance[final_report_sp_Inheritance['SZAreportCategory'].isin(['14','15','9','10'])]

#Those line is to only keep selected columns based on the apps need
#inal_report_general_2 = final_report_general_2[["Target.group","Disease","ClinVar_CLNDN","Genes","SZAID","Existing_variation","VARIANT_CLASS","Consequence","ClinVar_CLNHGVS","#CHROM","POS","ClinVar","REF","ALT","ClinVar_CLNSIG","Inheritances","DiseaseInfo","Zygosity","ClinVar_GENEINFO"]]
#final_report_general_3 = final_report_sp_Inheritance[final_report_sp_Inheritance['SZAreportCategory'].isin(['6','7','8','9','10','11','12','13','14','15'])]
#final_report_general_3 = final_report_general_3[["Target.group","Disease","ClinVar_CLNDN","Genes","SZAID","Existing_variation","VARIANT_CLASS","Consequence","ClinVar_CLNHGVS","#CHROM","POS","ClinVar","REF","ALT","ClinVar_CLNSIG","Inheritances","DiseaseInfo","Zygosity","ClinVar_GENEINFO"]]
final_report_general_4 = final_report_sp_Inheritance

#final_report_general_1.to_csv(outFile.replace('.txt','_sp_Inheritance_1.txt'), sep = "\t",index = None, header=True)
#final_report_general_2.to_csv(outFile.replace('.txt','_sp_Inheritance_2.txt'), sep = "\t",index = None, header=True)
#final_report_general_3.to_csv(outFile.replace('.txt','_sp_Inheritance_3.txt'), sep = "\t",index = None, header=True)
final_report_general_4.to_csv(outFile.replace('.txt','_sp_Inheritance_4.txt'), sep = "\t",index = None, header=True)

print("the file named _sp_Inheritance generated: XXX_sp_Inheritance_1.txt,XXX_sp_Inheritance_2.txt,XXX_sp_Inheritance_3.txt,XXX_sp_Inheritance_4")


#%% Cell 11 def append_data function
colNames = list(final_report_sp_Inheritance.columns)
def append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID,iinh, imark, imim, iExisting_variation, iCHROM, iPOS,
                                          iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                          iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory,istar,icompound_het,iSymbol):
    line_list.append(iline)
    TargetGroup.append(iTargetGroup)
    Disease.append(iDisease)
    Gene.append(iGene)
    SZAID.append(iSZAID)
    inh.append(iinh)
    mark.append(imark)
    mim.append(imim)
    Existing_variation.append(iExisting_variation)
    CHROM.append(iCHROM)
    POS.append(iPOS)
    Zygosity.append(iZygosity)
    REF.append(iREF)
    ALT.append(iALT)
    ID.append(iID)
    ClinVar_CLNSIG.append(iClinvar_CLNSIG)
    VARIANT_CLASS.append(iVARIANT_CLASS)
    Consequence.append(iConsequence)
    ClinVar_CLNHGVS.append(iClinVar_CLNHGVS)
    Inheritances.append(iheritance)
    VariantINFO.append(iVariantINFO)
    DiseaseInfo.append(iDiseaseInfo)
    SZAreportCategory.append(iSZAreportCategory)
    ReviewStar.append(istar)
    Compound_heterozygous.append(icompound_het)
    Symbol.append(iSymbol)


# create lists
TargetGroup = []
Disease = []
Gene = []
Symbol = []
SZAID = []
inh = []
mark = []
mim = []
Existing_variation =[]
VARIANT_CLASS =[]
Consequence =[]
ClinVar_CLNHGVS=[]
#Genotype = []
CHROM = []
POS = []
ID = []
REF=[]
ALT=[]
Zygosity = []
ConfidenceLevel = []
#SZAdiseaseID = []
ClinVar_CLNSIG=[]
Inheritances =[]
VariantINFO =[]
DiseaseInfo =[]
line_list =[]
VariantINFO =[]
SZAreportCategory =[]
ReviewStar = []
Compound_heterozygous = []


#%%Cell 12 def Remove the duplicated features fuction
def rm_duplicates(dupfile,nodupfile):
    with open(dupfile,'r') as f:
        #skip the title
        next(f)
        for line in f:
            line = line.replace('\n', '')  
            temp = line.split('\t') 
            iline = line 

            iClinVar_GENEINFO = temp[colNames.index('ClinVar_GENEINFO')]
            iClinVar_GENEINFO = [s.split(':')[0] for s in iClinVar_GENEINFO.split("&")]
            if len(iClinVar_GENEINFO) > 1:
                iClinVar_GENEINFO = list(set(iClinVar_GENEINFO))
                iClinVar_GENEINFO = ";".join(iClinVar_GENEINFO)
            iSZAID = temp[colNames.index('SZAID')]
            iSymbol = temp[colNames.index('SYMBOL')]
            iGene = temp[colNames.index('Genes')]
            iDiseaseInfo = temp[colNames.index('DiseaseInfo')]
            iinh = temp[colNames.index('Inheritance')]
            imark = temp[colNames.index('Mark')]
            #imim = temp[colNames.index('MIM')]
            imim = temp[colNames.index('MIM')].split(".")[0]
            #0807:
            iSZAreportCategory = temp[colNames.index('SZAreportCategory')]
            istar = temp[colNames.index('ClinVar_ReviewStar')]
            #new
            iREF=temp[colNames.index('REF')]
            iALT=temp[colNames.index('ALT')]
            iVARIANT_CLASS=temp[colNames.index('VARIANT_CLASS')]
            iClinvar_CLNSIG=temp[colNames.index('ClinVar_CLNSIG')].split("&")[0]
            iID = temp[colNames.index('ID')]
            iheritance = temp[colNames.index('Inheritances')]
            if iGene == 'Intergenic':
                iDisease = temp[colNames.index('ClinVar_CLNDN')]
            else:
                iDisease = temp[colNames.index('Disease')]
            if iDisease == '':
                sys.exit()
                #continue
            iTargetGroup = temp[colNames.index('Target.group')]
            if iTargetGroup == '':
                iTargetGroup = 'Expanded_clinvar_mono_diseases'
        #           iGenotype = temp[colNames.index('Genotype')]
            iCHROM = temp[colNames.index('#CHROM')]
            iPOS = temp[colNames.index('POS')]
            iZygosity = temp[colNames.index('Zygosity')]
            icompound_het = temp[colNames.index('Compound_heterozygous')]
            iExisting_variation = temp[colNames.index('Existing_variation')]
            iVariantINFO = iCHROM + ":" + iPOS + " " + iREF + ">" + iALT
            iConsequence = temp[colNames.index('Consequence')]
            iClinVar_CLNHGVS = temp[colNames.index('ClinVar_CLNHGVS')]
            iVariantINFO = iCHROM + ":" + iPOS + " " + iREF + ">" + iALT
            iIndex1 = [i for i, x in enumerate(TargetGroup) if  iTargetGroup == x]
            iIndex2 = [i for i, x in enumerate(Disease) if iDisease == x]
            iIndex3 = [i for i, x in enumerate(Gene) if iGene == x]
            iIndex9 = [i for i, x in enumerate(Symbol) if iSymbol == x]
            iIndex4 = [i for i, x in enumerate(SZAID) if iSZAID == x]
            iIndex5 = [i for i, x in enumerate(SZAreportCategory) if iSZAreportCategory == x]
            iIndex6 = [i for i, x in enumerate(inh) if iinh == x]
            iIndex7 = [i for i, x in enumerate(mark) if imark == x]
            iIndex8 = [i for i, x in enumerate(mim) if imim == x]
            #iIndex10 = [i for i, x in enumerate(ReviewStar) if istar == x]
            iIndex12345678 = [x for x in iIndex1 if x in iIndex2 and x in iIndex3 and x in iIndex4 and x in iIndex5 and x in iIndex6 and x in iIndex7 and x in iIndex8 ]
            iIndex12456789 = [x for x in iIndex1 if x in iIndex2 and x in iIndex8 and x in iIndex4 and x in iIndex5 and x in iIndex6 and x in iIndex7 and x in iIndex9 ]


            if 'SZAvar' in iSZAID and (iGene in iClinVar_GENEINFO or iGene == iClinVar_GENEINFO) and ";" not in iGene:
            #if len(iIndex12345) == 0:
                if len(iIndex12345678) == 0  | len(iIndex12456789) == 0:
                #For all different panels, keep independent records
                    append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID, iinh, imark, imim, iExisting_variation, iCHROM, iPOS,
                                    iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                    iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory, istar,icompound_het, iSymbol)

            elif 'SZAvar' in iSZAID  and ";" in iGene:
                if isinstance(iGene, list):
                    iGene = ';'.join(iGene)
                genes_set = set(iGene.split(';'))
                
                if isinstance(iClinVar_GENEINFO, list):
                    iClinVar_GENEINFO = ';'.join(iClinVar_GENEINFO)
                geneinfo_set = set(iClinVar_GENEINFO.split(';'))
                        
                if genes_set.issubset(geneinfo_set) or geneinfo_set.issubset(genes_set):
                    if len(iIndex12345678) == 0  | len(iIndex12456789) == 0:
                        append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID, iinh, imark, imim, iExisting_variation, iCHROM, iPOS,
                                    iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                    iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory, istar,icompound_het, iSymbol)

            elif 'Novel' in iSZAID:
                if len(iIndex12345678) == 0  | len(iIndex12456789) == 0 :
                    append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID,iinh, imark, imim,  iExisting_variation, iCHROM, iPOS,
                                    iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                    iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory, istar,icompound_het, iSymbol)
            
            elif "HGMD" in iClinvar_CLNSIG:
                if len(iIndex12345678) == 0  | len(iIndex12456789) == 0 :
                    append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID,iinh, imark, imim,  iExisting_variation, iCHROM, iPOS,
                                    iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                    iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory, istar,icompound_het,iSymbol)

            

#write nodup file
    with open(nodupfile,'w') as f:
        #print('okkk')
        colNames_new = ["final_target_group"] + colNames
        f.write("\t".join(colNames_new)+"\n")
        #f.write(line + "\n")
        #f.write('\t'.join(['Target group','Sub-group','Phenotype','Gene','SZAvarID','RS ID','VariantType','VariantTrans','HGVS','VariantINFO','CHROM','POS','ID','REF','ALT','NOTE','Disease information','Prevalence','Inheritance'])+'\n')
        for j in range(0,len(Disease)):
            if TargetGroup[j] in ['Health predipositions/Disease risk','Carrier-screening','Heriditary-cancer risk syndrome','Newborn-screening'] :
                #f.write('Basic (for healthy subjects)\tUsually used for:'+ TargetGroup[j]+'\t'+Disease[j]+'\t'+Gene[j]+'\t'+SZAID[j]+'\t'+ Existing_variation[j] +'\t' +VARIANT_CLASS[j]+ '\t' + Zygosity[j] +'\t'+ ClinVar_CLNHGVS[j]+ '\t' + VariantINFO[j] +'\t' +CHROM[j]+'\t'+POS[j]+'\t'+ ID[j] +'\t'+ REF[j]+'\t'+ ALT[j]+'\t'+ ClinVar_CLNSIG[j]+'\t'+ DiseaseInfo[j]+'\t'+ 'Not_provided'+ '\t'+ Inheritances[j] +'\n')
                f.write('Basic (for healthy subjects),Usually used for:'+ TargetGroup[j]+'\t' + line_list[j] +"\n")
                #f.write('Basic (for healthy subjects)\tUsually used for:'+TargetGroup[j]+'\t'+Disease[j]+'\t'+Gene[j]+'\t'+SZAID[j]+'\t'+'Not provided\tNot provided\tNot provided\tNot provided\t'+Genotype[j]+'\t'+CHROM[j]+'\t'+POS[j]+'\t'+'Not provided\tNot provided\tNot provided\t'+'Pathogenic/Likely pathogenic'+'\tTBD\tTBD\n')
            elif TargetGroup[j] in ['Expanded_clinvar_mono_diseases','Expanded_mono_rare_diseases']:
                #f.write('Extended (for potential patients)\t'+ TargetGroup[j]+'\t'+Disease[j]+'\t'+Gene[j]+'\t'+SZAID[j]+'\t'+ Existing_variation[j] +'\t' +VARIANT_CLASS[j]+ '\t' + Zygosity[j] +'\t'+ ClinVar_CLNHGVS[j]+ '\t' + VariantINFO[j] +'\t' +CHROM[j]+'\t'+POS[j]+'\t'+ ID[j] +'\t'+ REF[j]+'\t'+ ALT[j]+'\t'+ ClinVar_CLNSIG[j]+'\t'+ DiseaseInfo[j]+'\t'+ 'Not_provided'+ '\t'+ Inheritances[j] +'\n')
                f.write('Extended (for potential patients),' + TargetGroup[j]+'\t' + line_list[j] + "\n")

#%%cell 13 Generate output files 

# input files

dupfile = outFile.replace('.txt','_sp_Inheritance_4.txt')
nodupfile = outFile.replace('.txt','_sp_Inheritance_4_nodup.txt')

rm_duplicates(dupfile,nodupfile) # The function directly remove duplicated and write the output file
print("generated no dup file,", nodupfile)

#dupfile2 = outFile.replace('.txt','_sp_Inheritance_3.txt')
#nodupfile2 = outFile.replace('.txt','_sp_Inheritance_3_nodup.txt')
#rm_duplicates(dupfile2,nodupfile2)
#print("generated no dup file,", nodupfile2)



# if Upload_To_Cloud == True:

#     ssh = paramiko.SSHClient()

#     ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

#     ssh.connect('ec2-16-171-6-238.eu-north-1.compute.amazonaws.com', username='ubuntu', password='AWuWqFRGHUhHV59!', port=22)

#     scp = SCPClient(ssh.get_transport())



#     remote_file_path = '/home/ubuntu/bckndsrv/RawDataOutputs(GenomicSZA)/' + nodupfile.split('/')[-1]



#     scp.put(nodupfile, remote_file_path)

#     print("File transferred successfully to AWS directory:", remote_file_path)



#     # Close the SCP and SSH sessions

#     scp.close()

#     ssh.close()


#%% Cell 14 HGMD variants merge
nodupfile_data = pd.read_csv(nodupfile, sep='\t')
nodupfile_data['POS'] = nodupfile_data['POS'].astype(str).str.split('.').str[0]
nodupfile_data['POS'] = nodupfile_data['POS'].apply(lambda x: int(float(x)) if pd.notna(x) else x)
nodupfile_data['ClinVar_ID'] = nodupfile_data['ClinVar_ID'].astype(str).str.split('.').str[0]
nodupfile_data['ClinVar_ReviewStar'] = nodupfile_data['ClinVar_ReviewStar'].astype(str).str.split('.').str[0]
nodupfile_data['MIM'] = nodupfile_data['MIM'].astype(str).str.split('.').str[0]
nodupfile_data['ClinVar'] = nodupfile_data['ClinVar'].astype(str).str.split('.').str[0]
nodupfile_data['Gene'] = nodupfile_data['Gene'].astype(str).str.split('.').str[0]
nodupfile_data['SZAreportCategory'] =nodupfile_data['SZAreportCategory'].apply(lambda x: int(float(x)) if pd.notna(x) else x)
nodupfile_data['POS'] = nodupfile_data['POS'].astype(str)

# Select relevant columns
col = HGMD_DM.columns[[0, 6, 7, 8, 9, 10, 11, 12, 13,14]]

# Create a matrix of NA values
a = pd.DataFrame(np.nan, index=np.arange(len(nodupfile_data)), columns=col)
nodupfile_data = pd.concat([nodupfile_data, a], axis=1)


for i in range(len(nodupfile_data)):

    ivar = nodupfile_data.loc[i, '#CHROM']+ '_' + nodupfile_data.loc[i, 'POS'] + '_' + nodupfile_data.loc[i, 'REF'] +'_' + nodupfile_data.loc[i, 'ALT']
    index = HGMD_DM[HGMD_DM['variant_info'] == ivar].index


    if len(index) > 0:
        df = HGMD_DM.loc[index]

        # changing data type 
        nodupfile_data.loc[i, 'gene'] = df['gene'].values[0]
        nodupfile_data.loc[i, "Reported.phenotype"] = df['Reported.phenotype'].values[0]
        nodupfile_data.loc[i, 'Variant.class'] = df['Variant.class'].values[0]
        nodupfile_data.loc[i, 'Mutation.Type'] = df['Mutation.Type'].values[0]
        nodupfile_data.loc[i, 'HGMD.accession'] = df['HGMD.accession'].values[0]
        nodupfile_data.loc[i, 'HGMD.codon.change'] = df['HGMD.codon.change'].values[0]
        nodupfile_data.loc[i, 'HGMD.amino.acid.change'] = df['HGMD.amino.acid.change'].values[0]
        nodupfile_data.loc[i, 'HGMD.splicing.mutation'] = df['HGMD.splicing.mutation'].values[0]
        nodupfile_data.loc[i, 'HGVS..nucleotide.'] = df['HGVS..nucleotide.'].values[0]
        nodupfile_data.loc[i, 'HGVS..protein.'] = df['HGVS..protein.'].values[0]


total_newline = pd.DataFrame()
for i in range(len(nodupfile_data)):
    hgmd_score = int(float(nodupfile_data.loc[i, "SZAreportCategory"]))
    dm_level = nodupfile_data.loc[i, "Variant.class"]
    
    if not pd.isna(nodupfile_data.loc[i,"Reported.phenotype"]):
        if dm_level == "DM":
            dm_score = 5 + hgmd_score
        elif dm_level == "DM?":
            dm_score = 4 + hgmd_score
        
        nodupfile_data.loc[i, "SZAreportCategory"] = dm_score


    # both above situation will run the following script
        hgmd_disease = nodupfile_data.loc[i,"Reported.phenotype"].split(" & ")
        nodup_disease = nodupfile_data.loc[i, "Disease"]
        nodup_disease = re.sub(r'[^\w\s]', ' ', nodup_disease)
        nodup_disease_word = set(nodup_disease.upper().split())
        
        all_newline = pd.DataFrame()
        count = 0
        for j in range(len(hgmd_disease)):
            disease_tmp = str(hgmd_disease[j])
            disease_tmp = re.sub(r'[^\w\s]', ' ', disease_tmp)
            disease_tmp_word = set(disease_tmp.upper().split())
            if disease_tmp_word.issubset(nodup_disease_word) or nodup_disease_word.issubset(disease_tmp_word) : # if diseases matched
                count += 1
                con_level = nodupfile_data.loc[i, "Gene.Disease.confidence.level"]

                if con_level == 'High_confidence':
                    final_score = 10 + dm_score
                elif con_level == 'Moderate_confidence':
                    final_score = 5 + dm_score
                elif con_level == 'Low_confidence':
                    final_score = dm_score # Do not exist
                newline = nodupfile_data.loc[i].copy()
                newline = pd.DataFrame(newline).transpose()
                newline['SZAreportCategory'] = final_score
                newline['Reported.phenotype'] = hgmd_disease[j]
                all_newline = pd.concat([all_newline, newline], axis=0, ignore_index=True)
                total_newline = pd.concat([total_newline, all_newline], axis=0, ignore_index=True)

            else: # if disease not matched 
                newline = nodupfile_data.loc[i].copy()
                newline = pd.DataFrame(newline).transpose()
                newline['SZAreportCategory'] = dm_score
                newline['Reported.phenotype'] = hgmd_disease[j]
                all_newline = pd.concat([all_newline, newline], axis=0, ignore_index=True)
                total_newline = pd.concat([total_newline, all_newline], axis=0, ignore_index=True)
        
        if len(all_newline) > 0 and count != 0: # if all matched, then replace the line with newline
            nodupfile_data = nodupfile_data.drop(index=i)

        
if len(total_newline) > 0 :
    nodupfile_data = pd.concat([nodupfile_data, total_newline], ignore_index=False)



cols_to_rename = [col for col in nodupfile_data.columns if 'P00' in col]
new_colsnames = {col: 'Sample_WGS' for col in cols_to_rename}
nodupfile_data.rename(columns=new_colsnames, inplace=True)


col_relocation = ['final_target_group', 'SZAID', 'Database', 'Database_type', 'SZAdiseaseID', 'ClinVar_CLNSIG', 'Variant.class','ClinVar_CLNDN',  'Reported.phenotype',
                  'SZAreportCategory', 'ClinVar_ReviewStar', 'ClinVar_ID', 'selCSQ', 'Zygosity', 'Genotype', 'Gene.Disease.confidence.level', 
                  'Target.group', 'Genes', 'Disease', 'Inheritance', 'Mark', 'MIM', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                    'FORMAT', 'Sample_WGS', 'Compound_heterozygous', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 
                    'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 
                    'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'REFSEQ_MATCH',
                'REFSEQ_OFFSET', 'GIVEN_REF', 'USED_REF', 'BAM_EDIT', 'SOURCE', 'SIFT', 'PolyPhen', 'HGVS_OFFSET', 'AF', 'gnomADe_AF', 
                'gnomADe_AFR_AF', 'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF',
                'gnomADe_SAS_AF', 'gnomADg_AF', 'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF', 'gnomADg_EAS_AF', 
                'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG',
                'SOMATIC', 'PHENO', 'ada_score', 'rf_score', 'REVEL', 'BayesDel_addAF_pred', 'BayesDel_addAF_rankscore', 'BayesDel_addAF_score',
                    'BayesDel_noAF_pred', 'BayesDel_noAF_rankscore', 'BayesDel_noAF_score', 'REVEL_rankscore', 'REVEL_score', 'SpliceAI_pred_DP_AG', 
                    'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 
                    'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL', 'am_class', 'am_genome', 'am_pathogenicity', 'am_protein_variant', 'am_transcript_id',
                    'am_uniprot_id', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info', 'ClinVar',
                    'ClinVar_CLNHGVS', 'ClinVar_CLNSIGINCL', 'ClinVar_CLNVC', 'ClinVar_GENEINFO', 'ClinVar_CLNDISDB', 'ClinVar_CLNSIGCONF', 'ClinVar_CLNREVSTAT',
                        'ClinVar_CLNDNINCL',  'Database_Type', 'Database_SZAID', 'Inheritances', 'DiseaseInfo', 'gene', 'Mutation.Type', 
                        'HGMD.accession', 'HGMD.codon.change', 'HGMD.amino.acid.change', 'HGVS..nucleotide.', 'HGVS..protein.']


nodupfile_data = nodupfile_data[col_relocation]
nodupfile_data.rename(columns = {"Variant.class": "HGMD.class", "Reported.phenotype": "HGMD.phenotype"}, inplace=True) 

nodupfile_data = nodupfile_data.drop_duplicates()
nodupfile_data.to_csv(nodupfile, sep="\t",index = None, header=True)
    

#%% Cell 15 GeneBe ACMG classification
import genebe as gnb

genebe_df = pd.read_csv(nodupfile, sep="\t",header=0)
#small_df = genebe_df.loc[:,["#CHROM","POS","REF","ALT"]]
#small_df = small_df.rename(columns={"#CHROM":"chr","POS":"pos","REF":"ref","ALT":"alt"})
#unique_small_df = small_df.drop_duplicates()
#annotated_df = gnb.annotate(unique_small_df,
#    genome='hg38',
#    use_ensembl=False,
#    use_refseq=True,
#    flatten_consequences=True,
#    output_format="dataframe")
#annotated_df = annotated_df.rename(columns={"chr":"#CHROM","pos":"POS","ref":"REF","alt":"ALT"})
#
#small_annotate_all = annotated_df[["#CHROM","POS","REF","ALT","gene_symbol",
#                                   "acmg_score","acmg_classification",'acmg_criteria']]
#small_annotate_all = small_annotate_all.rename(columns={'gene_symbol': 'Genes'})
#genebe_df = pd.merge(genebe_df, small_annotate_all, how="left", on=["#CHROM","POS","REF","ALT","Genes"])
#merged_df.to_csv(nodupfile, sep="\t",index=False)


#%% Cell 16 ClinGen gene-variants classification
clingene = pd.read_csv(clingen_file, sep="\t",header=0)
genebe_df["ClinVar_ID"] = genebe_df["ClinVar_ID"].str.split("&").str[0]
clingene = clingene[['#Variation', 'ClinVar Variation Id',
      'HGNC Gene Symbol', 'Disease', 'Mondo Id',
       'Mode of Inheritance', 'Assertion', 'Applied Evidence Codes (Met)',
       'Applied Evidence Codes (Not Met)', 'Summary of interpretation']]
clingene = clingene.rename(columns={"ClinVar Variation Id": "ClinVar_ID", "HGNC Gene Symbol":"Genes",
                                    "Assertion":"ClinGen_classification","Disease":"ClinGen_disease",
                                    "#Variation":"ClinGen_variant",
                                    "Applied Evidence Codes (Met)":"ClinGen_applied_evidence_codes"})
clingen_merge_df = pd.merge(genebe_df, clingene, how="left", on=["ClinVar_ID","Genes"])



#%%cell 17 Add GenCC and ClinGen gene-disease classification
genecc_clingen_classification = pd.read_csv(geneBaseFile, sep="\t",header=0)
genecc_clingen_classification.drop(columns='GenCC_classification_clingen', inplace=True)
clingen_merge_df = pd.merge(clingen_merge_df, genecc_clingen_classification, how="left", 
                            on=['Gene.Disease.confidence.level', 'Target.group', 'Genes', 'Disease',
                                'SZAdiseaseID', 'Inheritance', 'Mark', 'MIM', ])

#clingen_merge_df.to_csv(nodupfile, sep="\t",index = None, header=True) 


#%%cell 18 Add ontology  
ontology = pd.read_csv(ontology_file, sep=",",header = 0)
ontology = ontology[['Genes',"SZAdiseaseID","category"]]
clingen_merge_df = pd.merge(clingen_merge_df, ontology, on=['Genes','SZAdiseaseID'], how='left')
clingen_merge_df.to_csv(nodupfile, sep="\t",index = None, header=True) 

#%%cell 18 Simplify the columns

if run_simplify:


    s_colname = ['ClinVar_CLNSIG', 'HGMD.class',"acmg_classification",'acmg_criteria','ClinGen_classification','ClinGen_applied_evidence_codes','ClinVar_CLNDN', 'HGMD.phenotype',
                  'SZAreportCategory', 'ClinVar_ReviewStar', 'Zygosity', 'Genotype', 
                  'Gene.Disease.confidence.level', 'Target.group', 'Genes', 'Disease', 'Inheritance', 'Mark',
                  '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT', 'Sample_WGS', 'Compound_heterozygous',  'Consequence', 'IMPACT', 
                 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'Existing_variation',
                'SIFT', 'PolyPhen', 'HGVS_OFFSET', 'AF', 'gnomADe_AF', 'gnomADe_AFR_AF', 'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF', 
                 'gnomADe_SAS_AF', 'ada_score', 'rf_score', 'REVEL', 'BayesDel_addAF_pred', 
                'BayesDel_addAF_rankscore', 'BayesDel_addAF_score','BayesDel_noAF_pred', 'BayesDel_noAF_rankscore', 
                'BayesDel_noAF_score', 'REVEL_rankscore', 'REVEL_score', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 
                'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 
                    'SpliceAI_pred_DS_DL', 'am_class', 'am_genome', 'am_pathogenicity', 'am_protein_variant', 'am_transcript_id',
                    'am_uniprot_id', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info', 'Inheritances', 'DiseaseInfo', 'Mutation.Type', 
                        'HGMD.accession', 'HGMD.codon.change', 'HGMD.amino.acid.change', 'HGVS..nucleotide.', 'HGVS..protein.']

    simplified_nodupfile = nodupfile_data[s_colname]

# Only keep the rsID
    for i in range(len(simplified_nodupfile)):
        ids = str(simplified_nodupfile.loc[i, 'Existing_variation'])
        ids = ids.split('&')
        rsid = [id for id in ids if 'rs' in id]
        rsid = ''.join(rsid)
        simplified_nodupfile.loc[i, 'Existing_variation'] = rsid


    simplified_nodupfile['ClinVar_ReviewStar'] = simplified_nodupfile['ClinVar_ReviewStar'].astype(str).str.split('.').str[0]
    simplified_nodupfile.to_csv(outFile.replace('.txt','_simplified.txt'), sep = "\t",index = None, header=True)

