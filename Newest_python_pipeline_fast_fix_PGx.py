#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:41:09 2023

@author: xiyas
"""
# GenomeReportingPipeline version 2 on SZA sercer
# June 2023: Adding the male/female distinguishment for "Hemizygous"
# Jan 2024: vep gz file compatible
# April 2024: Traits report, if jText[colNames_CSQ.index('Existing_variation')] in trait_list: saveFlag5 ="Trait_var"
#   rs17822931&CM062373&COSV62323034 should be separated as rs17822931

#%%Cell 1 Define the path for parameters and databases
import sys
import pandas as pd
import numpy as np
import os
import gzip

# 1. choose which enviroment running local uppmax, SZA
path = 'local'

# 2. choose which cohort belongs to, Turkish or Swedish
#cohort = 'Turkish'
#print("Turkish vcf files")
# 3. choose whether need to run the GWAS and Pharmaco findings
run_GWAS = False
run_Pharmaco = False
run_traits = True

# Check the input value and execute code accordingly
if path == 'local':
    #fileName = sys.argv[1]
    #outFile = sys.argv[2]
    fileName = '/Users/xiyas/V2_Genome_reporting/Sample_data/P001_01.hard-filtered.vcf.gz_vep_annotated.vcf.gz'
    #fileName = '/Users/xiyas/Downloads/P0064_1_vep_annotated.vcf.gz'
    #fileName= '/Users/xiyas/V2_Genome_reporting/Sample_data/P001_10.hard-filtered.vcf.gz_vep_annotated.vcf'
    output_directory = "/Users/xiyas/V2_Genome_reporting/python_output_turkish_275/Other"
    outFile = os.path.join(output_directory, 'test_your_pipeline_traits.txt')

    # ##define the name for the output file
    geneBaseFile = '/Users/xiyas/V2_Genome_reporting/database-file-v2/GeneDB.txt'
    #diseaesDBfile = '/Users/xiyas/V2_Genome_reporting/database-file-v2/DiseaseDB_version_2_Beta_curated.txt'
    diseaesDBfile = '/Users/xiyas/V2_Genome_reporting/database-file-v2/diseaseDB_1115_3.txt'
    OMIM_inheritance_DBfile = '/Users/xiyas/V2_Genome_reporting/database-file-v2/pheno_OMIM_all.txt'
    variantSumFile = "/Users/xiyas/V2_Genome_reporting/database-file-v2/All_variant_summary2_ver2.txt"

    #Read disease and inheritance DB
    DiseaseDB = pd.read_csv(diseaesDBfile, sep="\t",encoding="ISO-8859-1")
    DiseaseDB= DiseaseDB.replace(np.nan,"No info")
    OMIM_Inheritance_DB = pd.read_csv(OMIM_inheritance_DBfile,sep="\t",dtype={"phenotypeMimNumber": str})
    OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")

    print("Set path for local")

    if run_GWAS == True:
        GWAS_dbfile = '/Users/xiyas/V2_Genome_reporting/database-file-v2/Merged_GWAS_vcf.txt'
        GWAS_db =pd.read_csv(GWAS_dbfile, sep="\t")
        GWAS_db= GWAS_db.replace(np.nan,"No info")
        print("GWAS matching will run")

    if run_Pharmaco == True:
        Pharma_dbfile = '/Users/xiyas/V2_Genome_reporting/database-file-v2/Merged_Pharma_vcf.txt'
        Pharma_db =pd.read_csv(Pharma_dbfile, sep="\t")
        Pharma_db= Pharma_db.replace(np.nan,"No info")
        print("Pharmaco matching will run")

    if run_traits == True:
        Trait_dbfile = '/Users/xiyas/V2_Genome_reporting/database-file-v2/Reports_genome_databases_traits_merged.txt'
        Trait_db =pd.read_csv(Trait_dbfile,sep="\t",encoding = "ISO-8859-1")
        Trait_db = Trait_db.dropna(axis=1, how='all')
        Trait_db= Trait_db.replace(np.nan,"No info")
        trait_list=Trait_db['variants'].to_list()
        print("Traits matching will run")


elif path == 'uppmax':
    # Code for alternative 2
    fileName = sys.argv[1]
    outFile = sys.argv[2]
    geneBaseFile = '/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/GeneDB.txt'
    diseaesDBfile = '/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/diseaseDB_1115_3.txt'
    OMIM_inheritance_DBfile = '/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/pheno_OMIM_all.txt'
    variantSumFile = '/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/All_variant_summary2_ver2.txt'

    #Read disease DB
    DiseaseDB = pd.read_csv(diseaesDBfile, sep="\t",encoding="ISO-8859-1")
    DiseaseDB= DiseaseDB.replace(np.nan,"No info")
    OMIM_Inheritance_DB = pd.read_csv(OMIM_inheritance_DBfile,sep="\t",dtype={"phenotypeMimNumber": str})
    OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")

    print("Set path for uppmax")

    if run_GWAS == True:
        GWAS_dbfile = '/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Merged_GWAS_vcf.txt'
        GWAS_db =pd.read_csv(GWAS_dbfile, sep="\t")
        GWAS_db= GWAS_db.replace(np.nan,"No info")
        print("GWAS matching will run")

    if run_Pharmaco == True:
        Pharma_dbfile = '/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Merged_Pharma_vcf.txt'
        Pharma_db =pd.read_csv(Pharma_dbfile, sep="\t")
        Pharma_db= Pharma_db.replace(np.nan,"No info")
        print("Pharmaco matching will run")

    if run_traits == True:
        Trait_dbfile = '/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Reports_genome_databases_traits_merged.txt'
        Trait_db =pd.read_csv(Trait_dbfile,sep="\t",encoding = "ISO-8859-1")
        Trait_db = Trait_db.dropna(axis=1, how='all')
        Trait_db= Trait_db.replace(np.nan,"No info")
        trait_list=Trait_db['variants'].to_list()
        print("Trait matching will run")


elif path == 'SZA':
    fileName = sys.argv[1]
    outFile = sys.argv[2]
    # ####This should be the old path for 144 server, for the 240 server it should be /mnt/SZAPORTAL/database-file-v2
    geneBaseFile = '/mnt/SZAPORTAL/database-file-v2/GeneDB.txt'
    diseaesDBfile = '/mnt/SZAPORTAL/database-file-v2/diseaseDB_1115_3.txt'
    OMIM_inheritance_DBfile = '/mnt/SZAPORTAL/database-file-v2/pheno_OMIM_all.txt'
    variantSumFile = '/mnt/SZAPORTAL/database-file-v2/All_variant_summary2_ver2.txt'

    DiseaseDB = pd.read_csv(diseaesDBfile, sep="\t",encoding = "ISO-8859-1")
    DiseaseDB= DiseaseDB.replace(np.nan,"No info")
    OMIM_Inheritance_DB = pd.read_csv(OMIM_inheritance_DBfile,sep="\t",dtype={"phenotypeMimNumber": str})
    OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")

    print("Set path for SZA")

    if run_GWAS == True:
        GWAS_dbfile = '/sza_data/GenomicPipeline/database-file-v2/Merged_GWAS_vcf.txt'
        GWAS_db =pd.read_csv(GWAS_dbfile, sep="\t")
        GWAS_db= GWAS_db.replace(np.nan,"No info")
        print("GWAS matching will run")

    if run_Pharmaco == True:
        Pharma_dbfile = '/sza_data/GenomicPipeline/database-file-v2/Merged_Pharma_vcf.txt'
        Pharma_db =pd.read_csv(Pharma_dbfile, sep="\t")
        Pharma_db= Pharma_db.replace(np.nan,"No info")
        print("Pharmaco matching will run")

    if run_traits == True:
        Trait_dbfile = '/sza_data/GenomicPipeline/database-file-v2/Reports_genome_databases_traits_merged.txt'
        Trait_db =pd.read_csv(Trait_dbfile,sep="\t",encoding = "ISO-8859-1")
        Trait_db = Trait_db.dropna(axis=1, how='all')
        Trait_db= Trait_db.replace(np.nan,"No info")
        trait_list=Trait_db['variants'].to_list()
        print("Traits matching will run")


##NEW update 2.07: I add a new annotation, Clinvar review status for varisnts. so the CSQ from 88 to 89 columns
colNames_CSQ =['Allele','Consequence','IMPACT','SYMBOL',
  'Gene','Feature_type','Feature','BIOTYPE','EXON',
  'INTRON','HGVSc','HGVSp','cDNA_position',
  'CDS_position','Protein_position','Amino_acids','Codons',
  'Existing_variation',
  'DISTANCE','STRAND','FLAGS','VARIANT_CLASS',
  'SYMBOL_SOURCE','HGNC_ID','CANONICAL','SOURCE',
  'SIFT','PolyPhen','HGVS_OFFSET',
  'AF','gnomADe_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_ASJ_AF',
  'gnomADe_EAS_AF','gnomADe_FIN_AF','gnomADe_NFE_AF','gnomADe_OTH_AF',
  'gnomADe_SAS_AF','gnomADg_AF','gnomADg_AFR_AF','gnomADg_AMI_AF',
  'gnomADg_AMR_AF','gnomADg_ASJ_AF','gnomADg_EAS_AF','gnomADg_FIN_AF',
  'gnomADg_MID_AF','gnomADg_NFE_AF','gnomADg_OTH_AF','gnomADg_SAS_AF',
  'MAX_AF','MAX_AF_POPS','CLIN_SIG','SOMATIC',
  'PHENO','ada_score','rf_score','REVEL',
  'BayesDel_addAF_pred','BayesDel_addAF_rankscore','BayesDel_addAF_score',
  'BayesDel_noAF_pred','BayesDel_noAF_rankscore','BayesDel_noAF_score',
  'REVEL_rankscore','REVEL_score',
  'SpliceAI_pred_DP_AG','SpliceAI_pred_DP_AL','SpliceAI_pred_DP_DG',
  'SpliceAI_pred_DP_DL','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL',
  'SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','SpliceAI_pred_SYMBOL',
  'ClinVar','ClinVar_ID','ClinVar_CLNSIG','ClinVar_CLNDN','ClinVar_CLNHGVS',
  'ClinVar_CLNSIGINCL','ClinVar_CLNVC','ClinVar_GENEINFO','ClinVar_CLNDISDB',
  'ClinVar_CLNSIGCONF','ClinVar_CLNREVSTAT',
  'Database','Database_Type','Database_SZAID']

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
    check_male_flag = 0>1
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
        saveFlag1, saveFlag2,saveFlag3,saveFlag4, saveFlag5 = False, False, False, False, False

        ##this is going to analyse all transcript result and any one satisfied the criteria, the whole variant goes into report A/B.
        ## 1.26 2024: In the conference I checked rs6025, which should be SZAvar360812 a famous PGx genes for F5 gene. However it only appears on nodup4files, as a clinVar genes,
        ## not in PGx reports. so I realized this elif here should be all changed to "if", because elif ignores other cases if the first one satisfied the criteria
        for j in range(0,len(iText)):
            jText = iText[j].split('|')
            if jText[colNames_CSQ.index('Database_SZAID')] != '' and 'ClinP_LP_var' in jText[colNames_CSQ.index('Database_Type')].split('&'):
                #sys.exit('!')
                saveFlag1 = 'ClinP_LP_var'
            #if not saveFlag1 and ((jText[colNames_CSQ.index('MAX_AF')] == '') or (float(jText[colNames_CSQ.index('MAX_AF')]) < 0.05)):
            #2024.0311 update: gnomad Version 4.0 allele frequencies by VEP annotation =========if is old annotaion with gnomad 3, don't run this
            ## step 1 : get the MAX_AF --------------------------------------
            #start_index = colNames_CSQ.index('gnomAD_4_AF')
            #end_index = colNames_CSQ.index('gnomAD_4_AF_remaining') + 1
            #if jText[colNames_CSQ.index('MAX_AF')] == '' :
            #    jText[colNames_CSQ.index('MAX_AF')] = 0
            #gnomad_values = [float(val) if val != '' else 0 for val in jText[start_index:end_index]]
            #max_gnomad_value = max(gnomad_values)
            #max_af_value = jText[colNames_CSQ.index('MAX_AF')]
            #max_value = max(float(max_gnomad_value), float(max_af_value))
            #print(max_af_value)
            #print(max_gnomad_value)
            ## step 2: use the new MAX_AF from the comparison
            # ------------------------------------------------------------
            if not saveFlag1 and (jText[colNames_CSQ.index('MAX_AF')] == '' or float(jText[colNames_CSQ.index('MAX_AF')])<0.05):
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
                        jText[colNames_CSQ.index('BayesDel_noAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_noAF_score')])>-0.0570105)):
                    saveFlag2 ="Putative_var"
                #MAX_AF to gnomad global AF, and new filters need to be applied here too
            if jText[colNames_CSQ.index('Database_SZAID')] != '' and 'GWAS_var' in jText[colNames_CSQ.index('Database_Type')].split('&'):
                saveFlag3 = "GWAS_var"
            if jText[colNames_CSQ.index('Database_SZAID')] != '' and 'Pharma_var' in jText[colNames_CSQ.index('Database_Type')].split('&'):
                saveFlag4 = "Pharma_var"
            if run_traits == True:
                #################Bug 0401 ; Existing_variation can be rs174547&COSV63520556, need to split by & and then match the rsid db.
                if jText[colNames_CSQ.index('Existing_variation')].split('&')[0] in trait_list:
                    saveFlag5 ="Trait_var"
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
    if "chrY" in tLine:
            #print("Sex:Male")
        check_male_flag ="Male"
#    if "27771923" in tLine:
#        sys.exit("!")
    if i%100000 == 0:
        print(str(i)+' lines processed!')
    tLine = file.readline()

gender = "Male" if check_male_flag else "Female"

if check_male_flag :
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
with open(geneBaseFile,'r') as f:
    for line in f:
        line = line.replace('\n','')
        temp = line.split('\t')
        ConfidenceLevel.append(temp[0])
        TargetGroup.append(temp[1])
        geneBasedRef.append(temp[2])
        Disease.append(temp[3])
        SZAdiseaseID_GDB.append(temp[4])
geneBasedRefHeadings = [ConfidenceLevel[0],TargetGroup[0],geneBasedRef[0],Disease[0]]
geneBasedRef.pop(0)
Disease.pop(0)
ConfidenceLevel.pop(0)
TargetGroup.pop(0)
SZAdiseaseID_GDB.pop(0)

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


newContentHeadings = ['SZAID','Database','Database_type','SZAdiseaseID','ClinVar_CLNSIG','ClinVar_CLNDN','SZAreportCategory','ClinVar_ID','selCSQ','Zygosity','Genotype']+geneBasedRefHeadings+headings;
newContent = '\t'.join(newContentHeadings)+'\n'


#%%Cell 4 report A is for Clinvar Variants
# Part A

for i in range(0,len(reportA)):
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
        else:
            iZygo = 'Compound heterozygous'
        if check_male_flag:
            if "chrX" in reportA[i] and iZygo == 'Heterozygous':
                iZygo = "Hemizygous"
        for j in range(0,len(iText)):
            jText = iText[j].split('|')
            jGenes1 = jText[3]
            jGenes2Temp = jText[41].split('&')
            jGenes2 = []
            for k in range(0,len(jGenes2Temp)):
                kTemp = jGenes2Temp[k].split(':')
                jGenes2.append(kTemp[0])
            jGenes = list(set([jGenes1]+jGenes2))
    ######### To clear GWAS and pharmaco ######
    ###The reason to do this: GWAS,pharmaco, Clinvar can be overlapped, so several SZAvar can be actually one, just showing in both clinvar and gwas or ..
    ###This will have issues like ValueError: 'ClinP_LP_var' is not in list
    ###Is because in one of the transcripts, the way of express alternative allele is different with what clinvar used.
    ## expression of "-" will annotated Clinvar_P_LP, while "TGCTGC"
            #print(jText[colNames_CSQ.index('Database_Type')].split('&'))
            if jText[colNames_CSQ.index('Database_Type')]!= '':
                jCLINVARind = jText[colNames_CSQ.index('Database_Type')].split('&').index('ClinP_LP_var')
                jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jCLINVARind]
                #jDatabase = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
                jCLNID = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
                jSZAvarID = jText[colNames_CSQ.index('Database_SZAID')].split('&')[jCLINVARind]
                #if jSZAvarID = 'SZAvar178179':
        ######### To clear GWAS and pharmaco #########
                jCLNSIG = jText[colNames_CSQ.index('ClinVar_CLNSIG')]
                #print(jCLNSIG)
                #jCLNSIG = jText[colNames_CSQ.index('ClinVar_CLNSIG')].split('&')
                #print(jCLNSIG)
                jClinVar_CLNDN =  jText[colNames_CSQ.index('ClinVar_CLNDN')].replace('&_',',_').replace('&','|').replace('_',' ')
                if 'Pathogenic' in jCLNSIG:
                    scoreFlag = 5
                elif 'Likely_pathogenic' in jCLNSIG:
                    scoreFlag = 4
                #elif 'risk_factor' in jCLNSIG:
                #    scoreFlag = 1
                elif 'Conflicting_interpretations_of_pathogenicity' in jCLNSIG:
                    scoreFlag = 1
                else:
                    scoreFlag = 0
                #if 'KRT75' in jGenes:
                #    sys.exit('!')
                if scoreFlag > 0 and any([x for x in jGenes if x in geneBasedRef]): # if the associated gene can be mapped to GeneDB

                    jGeneList = [x for x in jGenes if x in geneBasedRef]
                    jInds = []
                    #if 'HBG2' in jGeneList:
                    #    sys.exit('!')
                    for k in range(0,len(jGeneList)):
                        kTemp = jGeneList[k]
                        jInds = jInds +[ii for ii, x in enumerate(geneBasedRef) if x == kTemp]

                    ######### test
                    kClinVarDiseases = jClinVar_CLNDN.split('|')
                    #print(kClinVarDiseases)
                    jVarSZAdiseaseIDs = []
                    jVarDiseaseNames = []
                    for k in range(0,len(kClinVarDiseases)):
                        kClinVarDisease_temp = kClinVarDiseases[k]
                        if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1:
                            kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]
                            jVarSZAdiseaseIDs.append(DiseaseID_DSDB[kIndex[0]])
                            jVarDiseaseNames.append(DiseaseName_DSDB[kIndex[0]])
                    #sys.exit('!')
                    ## problem here!!! (should match Clinvar disease type?)
                    ## fixed by if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs: SZAdiseaseID_GDB is the disease from GeneDB, to see
                    # whether it's matching with the ClinVar_CLNDN disease jVarSZAdiseaseIDs
                    for k in range(0,len(jInds)):
                        if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs:
                            #SZAscore =0
                            #print(SZAscore)
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
                            newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],Disease[jInds[k]]]+iTemp)+'\n'
                        else:
                            ########here, if SZAdiseaseID_GDB[jInds[k]] in jVarSZAdiseaseIDs is false, then no SZAscore value
                            SZAscore = scoreFlag
                            jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                            newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],geneBasedRef[jInds[k]],Disease[jInds[k]]]+iTemp)+'\n'

                elif scoreFlag >0 and any(jGenes):
                    #print("here")
                    #sys.exit('Found a gene set with '+', '.join(jGenes)+' that is not included in GeneDB ')
                    SZAscore = scoreFlag
                    kClinVarDiseases = jClinVar_CLNDN.split('|')
                    for k in range(0,len(kClinVarDiseases)):
                        kClinVarDisease_temp = kClinVarDiseases[k]
                        if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1:
                            kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]
                            jSZADiseaseID = DiseaseID_DSDB[kIndex[0]]
                            jDisease = DiseaseName_DSDB[kIndex[0]]
                        else:
                            jSZADiseaseID = ''
                            jDisease = jClinVar_CLNDN
                            #sys.exit('!')
                    newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases',';'.join(jGenes),jDisease]+iTemp)+'\n'
                elif scoreFlag >0:
                    SZAscore = scoreFlag
                    kClinVarDiseases = jClinVar_CLNDN.split('|')
                    for k in range(0,len(kClinVarDiseases)):
                        kClinVarDisease_temp = kClinVarDiseases[k]
                        if len([ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]) == 1:
                            kIndex = [ii for ii, x in enumerate(DiseaseName_DSDB) if kClinVarDisease_temp ==  x.replace(',','')]
                            jSZADiseaseID = DiseaseID_DSDB[kIndex[0]]
                            jDisease = DiseaseName_DSDB[kIndex[0]]
                        else:
                            jSZADiseaseID = ''
                            jDisease = jClinVar_CLNDN
                            #sys.exit('!')
                    #sys.exit('!')
                    newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jSZADiseaseID,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,'Low_confidence','Expanded_clinvar_mono_diseases','Intergenic',jDisease]+iTemp)+'\n'
                    #newContent = newContent+'\t'.join([jSZAvarID,jCLNID,jDatabaseType,jCLNSIG,jClinVar_CLNDN,str(SZAscore),jCLNID,iText[j],iZygo,iGenotype,'','','Intergenic','']+iTemp)+'\n'


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


#%% Cell 6 report D for pharmaco-annotation

if run_Pharmaco == True:

    from datetime import datetime
    now = datetime.now()
    newLineHeadings = ['Drug.s.','Genotype','Zygosity','pharma_genotype','SZAvarID','Annotation.Text','Score','Phenotype.Category','Phenotype.s.','Gene','Level.of.Evidence','Level.Modifiers','URL'];
    newLine = '\t'.join(newLineHeadings)+'\n'
    with open(outFile.replace('.txt','_pharmaco.txt'),'w') as f:
        f.write(newLine)
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportD)) + ' lines of pharmaco variants' )
    # Pharma part
    for i in range(0,len(reportD)):
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
                jPharmaind = jText[colNames_CSQ.index('Database_Type')].split('&').index('Pharma_var')
                jDatabaseType = jText[colNames_CSQ.index('Database_Type')].split('&')[jPharmaind]
                #jDatabase = jText[colNames_CSQ.index('Database')].split('&')[jCLINVARind]
                jPharmaID = jText[colNames_CSQ.index('Database')].split('&')[jPharmaind]
                jSZAvarID = jText[colNames_CSQ.index('Database_SZAID')].split('&')[jPharmaind]
                ######### To match the Pharma_DB
                ### The pharma_genotype for one SZAID is multiple!!
                # so need to use genotype_ind
                jpharma_genotype_list = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['pharma_genotype'].to_list()
                if iGenotype in jpharma_genotype_list:
                    genotype_ind = jpharma_genotype_list.index(iGenotype)
                    jpharma_genotype = jpharma_genotype_list[genotype_ind]
                    jDrug= Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Drug.s.'].to_list()[genotype_ind]
                    jAnnotation = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Annotation.Text'].to_list()[genotype_ind]
                    jScore = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Score'].to_list()[genotype_ind]
                    jPheno_category = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Phenotype.Category'].to_list()[genotype_ind]
                    jPheno = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Phenotype.s.'].to_list()[genotype_ind]
                    jGene = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Gene'].to_list()[genotype_ind]
                    jevidence = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Level.of.Evidence'].to_list()[genotype_ind]
                    jlevel_modifier = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['Level.Modifiers'].to_list()[genotype_ind]
                    jURL = Pharma_db[Pharma_db['SZAID']== jSZAvarID]['URL'].to_list()[genotype_ind]
                    newLine_temp ='\t'.join([jDrug,iGenotype,iZygo,jpharma_genotype,jSZAvarID,jAnnotation,str(jScore),jPheno_category,jPheno,jGene,jevidence,jlevel_modifier,jURL])+'\n'
                    with open(outFile.replace('.txt','_pharmaco.txt'),'a') as f:
                        f.write(newLine_temp)
        #newLineHeadings = ['Drug(s)','Genotype','pharma_genotype','SZAvarID','Annotation Text','Score','Phenotype Category','Phenotype(s)','Gene','Level of Evidence','Level Modifiers','URL'];

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
        if check_male_flag:
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
                if (jText[colNames_CSQ.index('IMPACT')] == 'HIGH' or (
                        jText[colNames_CSQ.index('ada_score')] != '' and float(jText[colNames_CSQ.index('ada_score')])>0.6) or (
                        jText[colNames_CSQ.index('rf_score')] != '' and float(jText[colNames_CSQ.index('rf_score')])>0.6) or (
                        jText[colNames_CSQ.index('REVEL')]!= '' and float(jText[colNames_CSQ.index('REVEL')])>0.75) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DG')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_DL')])>0.5) or (
                        jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')]!= '' and float(jText[colNames_CSQ.index('SpliceAI_pred_DS_AG')])>0.5) or (
                        jText[colNames_CSQ.index('BayesDel_addAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_addAF_score')])>0.0692655) or (
                        jText[colNames_CSQ.index('BayesDel_noAF_score')]!= '' and float(jText[colNames_CSQ.index('BayesDel_noAF_score')])>-0.0570105)):

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
                                sys.exit('Unknown confidence score level for the report!')
                            jSZADiseaseID = SZAdiseaseID_GDB[jInds[k]]
                            #print(jGenes)
                            #SZAscore = 2+('High_confidence' == ConfidenceLevel[jInds[k]])*5+(jText[colNames_CSQ.index('IMPACT')]=='HIGH')*1
                            newContent = newContent+'\t'.join(
                                ['NovelTrans'+headings[9]+'_'+str(varCount+1),'NA','NA',
                                  jSZADiseaseID,'NA','NA',str(SZAscore),'NA',iText[j],iZygo,
                                  iGenotype,ConfidenceLevel[jInds[k]],TargetGroup[jInds[k]],
                                  geneBasedRef[jInds[k]],Disease[jInds[k]]]+iTemp)+'\n'
    if varFlag == 1:
        varCount = varCount+1

print("putative variants matched the GeneDB:",varCount)

# write final report
outfid = open(outFile,'w')
outfid.write(newContent)
outfid.close()

print("outFile generated ")

#%% Cell 9 length file generation
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

# print(f'Report lengths saved to {length_filename}')
# print(f"Input File: {fileName}")
# print(f"Gender: {gender}")
# print(f"Lengths saved to '{length_filename}'")

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
final_report_sp[colNames_CSQ]=final_report_sp['selCSQ'].apply(lambda x: pd.Series(str(x).split("|")))
# save as new format
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
        iGene = temp[13]
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
print("now only for file/4")
###Remove unneeded file
###point : there are some varints didn't match the SZAdiseaseID. temp[3] is empty.
#os.remove(outFile_sp_Inheritance) # As same as final_report_general_4,_sp_Inheritance_4.txt
#print(f"{outFile_sp_Inheritance} has been successfully deleted.")

#%% Cell 11 def append_data function

colNames = list(final_report_sp_Inheritance.columns)
def append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID, iExisting_variation, iCHROM, iPOS,
                                          iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                          iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory):
    line_list.append(iline)
    TargetGroup.append(iTargetGroup)
    Disease.append(iDisease)
    Gene.append(iGene)
    SZAID.append(iSZAID)
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


# create lists
TargetGroup = []
Disease = []
Gene = []
SZAID = []
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


#%%Cell 12 def Remove the duplicated features fuction
def rm_duplicates(dupfile,nodupfile):
    with open(dupfile,'r') as f:
        #skip the title
        next(f)
        for line in f:
            line = line.replace('\n','')
            temp = line.split('\t')
            iline = line
            iClinVar_GENEINFO = temp[colNames.index('ClinVar_GENEINFO')]
            iClinVar_GENEINFO = [s.split(':')[0] for s in iClinVar_GENEINFO.split("&")]
            #iClinVar_GENEINFO = iClinVar_GENEINFO.split(":")[0]

            iSZAID = temp[colNames.index('SZAID')]
            iGene = temp[colNames.index('Genes')]
            iDiseaseInfo = temp[colNames.index('DiseaseInfo')]
            #0807:
            iSZAreportCategory = temp[colNames.index('SZAreportCategory')]
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
            iExisting_variation = temp[colNames.index('Existing_variation')]
            iVariantINFO = iCHROM + ":" + iPOS + " " + iREF + ">" + iALT
  #          iConfidenceLevel = temp[colNames.index('SZAreportCategory')]
  #          iSZAdiseaseID = temp[colNames.index('SZAdiseaseID')]
            iConsequence = temp[colNames.index('Consequence')]
            iClinVar_CLNHGVS = temp[colNames.index('ClinVar_CLNHGVS')]
            iVariantINFO = iCHROM + ":" + iPOS + " " + iREF + ">" + iALT
            iIndex1 = [i for i, x in enumerate(TargetGroup) if  iTargetGroup == x]
            iIndex2 = [i for i, x in enumerate(Disease) if iDisease == x]
            iIndex3 = [i for i, x in enumerate(Gene) if iGene == x]
            iIndex4 = [i for i, x in enumerate(SZAID) if iSZAID == x]
            iIndex5 = [i for i, x in enumerate(SZAreportCategory) if iSZAreportCategory == x]

            iIndex12345 = [x for x in iIndex1 if x in iIndex2 and x in iIndex3 and x in iIndex4 and x in iIndex5]
            #iIndex2345 = [x for x in iIndex2 if x in iIndex3 and x in iIndex4 and x in iIndex5]
            #iIndex123 = [x for x in iIndex1 if x in iIndex2 and x in iIndex3]
            #check if Clinvar iGene in iClinvar_GENEINFO
            if 'SZAvar' in iSZAID and iGene in iClinVar_GENEINFO:
                if len(iIndex12345) == 0:
                #For all different panels, keep independent records
                    append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID, iExisting_variation, iCHROM, iPOS,
                                  iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                  iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory)

            elif 'Novel' in iSZAID:
                if len(iIndex12345) == 0:
                    append_data_to_lists(iline, iTargetGroup, iDisease, iGene, iSZAID, iExisting_variation, iCHROM, iPOS,
                                  iZygosity, iREF, iALT, iID, iClinvar_CLNSIG, iVARIANT_CLASS, iConsequence,
                                  iClinVar_CLNHGVS, iheritance, iVariantINFO, iDiseaseInfo, iSZAreportCategory)

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

# with open(outFile.replace('.txt','_sp_Inheritance_1_nodup_report_sza.txt'),'w') as f:
#     #colNames_new = ["final_target_group","Sub-group","VariantINFO"] + colNames #bunu ben yaptım sonra açarsın
#     #f.write(colNames + "\n")
#     #f.write("\t".join(colNames_new)+"\n") #bunuda ben ekledim
#     #f.write(line + "\n")
#     f.write('\t'.join(['Target group','Sub-group','Phenotype','Gene','SZAvarID','RS ID','VariantType','VariantTrans','HGVS','VariantINFO','CHROM','POS','ID','REF','ALT','NOTE','Disease information','Prevalence','Inheritance'])+'\n')
#     #f.write('\t'.join(['Target group','Sub-group','Phenotype','Gene','SZAvarID','RS ID','VariantType','VariantTrans','HGVS','VariantINFO','CHROM','POS','ID','REF','ALT','NOTE','Disease information','Prevalence','Inheritance'])+'\n')
#     for j in range(0,len(Disease)):
#         if TargetGroup[j] in ['Health predipositions/Disease risk','Carrier-screening','Heriditary-cancer risk syndrome','Newborn-screening'] and 'Novel' not in SZAID[j]:
#             f.write('Basic (for healthy subjects)\tUsually used for:'+ TargetGroup[j]+'\t'+Disease[j]+'\t'+Gene[j]+'\t'+SZAID[j]+'\t'+ Existing_variation[j] +'\t' +VARIANT_CLASS[j]+ '\t' + Zygosity[j] +'\t'+ ClinVar_CLNHGVS[j]+ '\t' + VariantINFO[j] +'\t' +CHROM[j]+'\t'+POS[j]+'\t'+ ID[j] +'\t'+ REF[j]+'\t'+ ALT[j]+'\t'+ ClinVar_CLNSIG[j]+'\t'+ DiseaseInfo[j]+'\t'+ 'Not_provided'+ '\t'+ Inheritances[j] +'\n')
#             #f.write('Basic (for healthy subjects)\tUsually used for:'+ TargetGroup[j]+'\t'+VariantINFO[j] +'\t' + line_list[j] + "\n")
#             #f.write('Basic (for healthy subjects)\tUsually used for:'+TargetGroup[j]+'\t'+Disease[j]+'\t'+Gene[j]+'\t'+SZAID[j]+'\t'+'Not provided\tNot provided\tNot provided\tNot provided\t'+Genotype[j]+'\t'+CHROM[j]+'\t'+POS[j]+'\t'+'Not provided\tNot provided\tNot provided\t'+'Pathogenic/Likely pathogenic'+'\tTBD\tTBD\n')
#         elif TargetGroup[j] in ['Expanded_clinvar_mono_diseases','Expanded_mono_rare_diseases'] and 'Novel' not in SZAID[j]:
#             f.write('Extended (for potential patients)\t'+ TargetGroup[j]+'\t'+Disease[j]+'\t'+Gene[j]+'\t'+SZAID[j]+'\t'+ Existing_variation[j] +'\t' +VARIANT_CLASS[j]+ '\t' + Zygosity[j] +'\t'+ ClinVar_CLNHGVS[j]+ '\t' + VariantINFO[j] +'\t' +CHROM[j]+'\t'+POS[j]+'\t'+ ID[j] +'\t'+ REF[j]+'\t'+ ALT[j]+'\t'+ ClinVar_CLNSIG[j]+'\t'+ DiseaseInfo[j]+'\t'+ 'Not_provided'+ '\t'+ Inheritances[j] +'\n')
#             #f.write('Extended (for potential patients)\t' + TargetGroup[j]+'\t' + VariantINFO[j] +'\t' +  line_list[j] + "\n"

#%%cell 13 main

# input files

dupfile = outFile.replace('.txt','_sp_Inheritance_4.txt')
nodupfile = outFile.replace('.txt','_sp_Inheritance_4_nodup.txt')

rm_duplicates(dupfile,nodupfile)
print("generated no dup file,", nodupfile)

#dupfile2 = outFile.replace('.txt','_sp_Inheritance_3.txt')
#nodupfile2 = outFile.replace('.txt','_sp_Inheritance_3_nodup.txt')
#rm_duplicates(dupfile2,nodupfile2)
#print("generated no dup file,", nodupfile2)

#%% Cell 14 Variant summary table

# TargetGroup_ref = []
# Disease_ref = []
# VarCount_ref = []
# Gene_ref = []
# with open(variantSumFile,'r') as f:
#     for line in f:
#         line = line.replace('\n','')
#         temp = line.split('\t')
#         TargetGroup_ref.append(temp[0])
#         Disease_ref.append(temp[1])
#         VarCount_ref.append(temp[2])
#         Gene_ref.append(temp[3])
# TargetGroup_ref.pop(0)
# Disease_ref.pop(0)
# VarCount_ref.pop(0)
# Gene_ref.pop(0)


# # generate variant summary table for the sample
# print("generating variant summary table for the sample")
# print("only the variants score as 9,10,14,15, Clinvar P/LP variants in HIGH/MODERATE confidence disases in the variant summary file")

# TargetGroup = []
# Disease = []
# Gene = []
# SZAvarID = []
# Genotype = []
# CHROM = []
# POS = []
# Zygosity = []
# ConfidenceLevel = []
# SZAdiseaseID = []
# count = 0
# tCount = 0
# VarCount_summary = [0]*len(VarCount_ref)

# #print(VarCount_summary)

# with open(outFile.replace('.txt', '_sp_Inheritance_3.txt'), 'r') as f:
#     for line in f:
#         line = line.replace('\n','')
#         temp = line.split('\t')
#         #iDiseaseClinVar = temp[4]
#         if temp[0] == 'SZAID':
#             colNames = temp
#         else:
#             iSZAID = temp[colNames.index('SZAID')]
#             iGene = temp[colNames.index('Genes')]
#             if iGene == 'Intergenic':
#                 iDisease = temp[colNames.index('ClinVar_CLNDN')]
#             else:
#                 iDisease = temp[colNames.index('Disease')]
#             if iDisease == '':
#                 sys.exit()
#             iTargetGroup = temp[colNames.index('Target.group')]
#             if iTargetGroup == '':
#                 iTargetGroup = 'Expanded_clinvar_mono_diseases'
#             iGenotype = temp[colNames.index('Genotype')]
#             iCHROM = temp[colNames.index('#CHROM')]
#             iPOS = temp[colNames.index('POS')]
#             iZygosity = temp[colNames.index('Zygosity')]
#             iConfidenceLevel = temp[colNames.index('SZAreportCategory')]
#             iSZAdiseaseID = temp[colNames.index('SZAdiseaseID')]
#             iClinVar_GENEINFO = temp[colNames.index('ClinVar_GENEINFO')]
#             iClinVar_GENEINFO = [s.split(':')[0] for s in iClinVar_GENEINFO.split("&")]
#             #iGENEINFO = [s.split(':')[0] for s in temp[colNames.index('selCSQ')].split('|')[-5].split('&')]
#             #iGENEINFO = [s.split(':')[0] for s in temp[colNames.index('selCSQ')].split('|')[-14].split('&')]
#             ######IF you change the column numbers of the selCSQ then you need change this line
#             #print(temp[colNames.index('selCSQ')].split('|')[-15].split('&'))
#             #print(temp[colNames.index('selCSQ')].split('|'))
#             #iGENEINFO = [s.split(':')[0] for s in temp[colNames.index('selCSQ')].split('|')[-15].split('&')]

#             #print(iGENEINFO)
#             iIndex1 = [i for i, x in enumerate(TargetGroup) if  iTargetGroup == x]
#             iIndex2 = [i for i, x in enumerate(Disease) if iDisease == x]
#             iIndex3 = [i for i, x in enumerate(Gene) if iGene == x]
#             iIndex4 = [i for i, x in enumerate(SZAvarID) if iSZAID == x]
#             iIndex124 = [x for x in iIndex1 if x in iIndex2 and x in iIndex4]
#             ###9,10,14,15: clinvar variants
#             # if iSZAID =="SZAvar9718":
#             #     sys.exit('!')
#             #if len(iIndex124) == 0 and iConfidenceLevel in ['9','10','14','15'] and iGene in iGENEINFO:
#             if len(iIndex124) == 0 and iConfidenceLevel in ['9','10','14','15'] and iGene in iClinVar_GENEINFO:
#                 #print('Clinvar_detected')
#                 jIndex1 = [ii for ii, xx in enumerate(TargetGroup_ref) if iTargetGroup == xx]
#                 jIndex2 = [ii for ii, xx in enumerate(Disease_ref) if iDisease == xx]
#                 jIndex12 = [xx for xx in jIndex1 if xx in jIndex2]
#                 #this is start of adding phenotypes detected by Clinvar
#                 TargetGroup.append(iTargetGroup)
#                 Disease.append(iDisease)
#                 Gene.append(iGene)
#                 SZAvarID.append(iSZAID)
#                 for j in jIndex12:
#                     #print(jIndex12)
#                     VarCount_summary[j] = VarCount_summary[j]+1
#                     #print(VarCount_summary[j])
#                     #print(iSZAID)


# with open(outFile.replace('.txt','_variantSummary1.txt'),'w') as f:
#     f.write('\t'.join(['Target group','Phenotype','Detected','Total','Gene','Novel'])+'\n')
#     for j in range(0,len(TargetGroup_ref)):
#         if TargetGroup_ref[j] in ['Health predipositions/Disease risk','Carrier-screening','Heriditary-cancer risk syndrome','Newborn-screening']:
#             f.write('Basic (for healthy subjects): usually used for '+TargetGroup_ref[j]+'\t'+Disease_ref[j]+'\t'+str(VarCount_summary[j])+'\t'+str(VarCount_ref[j])+'\t'+str(Gene_ref[j])+'\n')
# with open(outFile.replace('.txt','_variantSummary2.txt'),'w') as f:
#     f.write('\t'.join(['Target group','Phenotype','Detected','Total','Gene','Novel'])+'\n')
#     for j in range(0,len(TargetGroup_ref)):
#         if TargetGroup_ref[j] in ['Health predipositions/Disease risk','Carrier-screening','Heriditary-cancer risk syndrome','Newborn-screening']:
#             f.write('Basic (for healthy subjects): usually used for '+TargetGroup_ref[j]+'\t'+Disease_ref[j]+'\t'+str(VarCount_summary[j])+'\t'+str(VarCount_ref[j])+'\t'+str(Gene_ref[j])+'\n')
#         elif TargetGroup_ref[j] in ['Expanded_mono_rare_diseases','Expanded_clinvar_mono_diseases']:
#             f.write('Extended (for potential patients): '+TargetGroup_ref[j]+'\t'+Disease_ref[j]+'\t'+str(VarCount_summary[j])+'\t'+str(VarCount_ref[j])+'\t'+str(Gene_ref[j])+'\n')



# print("generated clinvar variant count files")
