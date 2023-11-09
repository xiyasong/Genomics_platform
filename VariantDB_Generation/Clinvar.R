library(qdapRegex)
#install.packages('R.utils')
library(R.utils)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)

#clin_vcf <- fread("/Users/xiyas/Report_genome/New_v2_database/clinvar_P_LP_RF_CLNSIG_1.vcf.gz")

#INFO_field <- strsplit(clin_vcf$INFO,';')


#GENEINFO <- lapply(INFO_field, grep, pattern = 'GENEINFO', value = TRUE)
#clin_vcf$GENEINFO <- GENEINFO
#clin_vcf$GENEINFO = gsub("GENEINFO=","",clin_vcf$GENEINFO)
#clin_vcf <- clin_vcf %>% separate(GENEINFO, c('gene_name', 'rm'),sep= ":")
#clin_vcf$rm <- NULL

#CLNSIG<- lapply(INFO_field, grep, pattern = 'CLNSIG', value = TRUE)
#clin_vcf$CLNSIG =  CLNSIG
#clin_vcf$CLNSIG = gsub("CLNSIG=","",clin_vcf$CLNSIG)

#RS <- lapply(INFO_field, grep, pattern = 'RS=', value = TRUE)
#clin_vcf$rs_id = RS
#clin_vcf$rs_id = gsub("RS=","",clin_vcf$rs_id)

#CLNDN <-lapply(INFO_field, grep, pattern = 'CLNDN', value = TRUE)
#clin_vcf$Phenotype = CLNDN
#clin_vcf$Phenotype = gsub("CLNDN=","",clin_vcf$Phenotype)

#clin_vcf_new <- clin_vcf
#clin_vcf_new$class <- rep("CLIN_P_LP_Rf")
#clin_vcf_new$QUAL <- NULL
#clin_vcf_new$FILTER <- NULL
#clin_vcf_new$id <- NULL
#clin_vcf_new$CLNSIG <- unlist(clin_vcf_new$CLNSIG)
#clin_vcf_new$rs_id <- unlist(clin_vcf_new$rs_id)
#clin_vcf_new$rs_id[clin_vcf_new$rs_id=="character(0)"] <- "NA"


#clin_vcf_new$rs_id <- paste0("rs",clin_vcf_new$rs_id)
#clin_vcf_new$INFO <- NULL

#------sza id clinvar file
library(dplyr)
clin_vcf_sza_id <- fread("/Users/xiyas/V2_Genome_reporting/clinvar_20221129.P.LP.modify.vcf.gz")

clin_vcf_sza_id <- clin_vcf_sza_id %>% distinct()

clin_vcf_sza_id$INFO <- rep("Type=ClinP_LP_var;")
clin_vcf_sza_id$INFO <- paste0(clin_vcf_sza_id$INFO,"SZAID=SZAvar")

#------------------Combine all the three database tables
colnames(clin_vcf_sza_id)[1] <- "CHROM"
Whole_database_GWAS_Phar_Clin <- rbind(clin_vcf_sza_id,Merged_GWAS_vcf_for_sza_id,Merged_pharma_vcf_for_sza_id)
Whole_database_GWAS_Phar_Clin$INFO <- paste0(Whole_database_GWAS_Phar_Clin$INFO, seq(1,361495))
dim(Whole_database_GWAS_Phar_Clin)
write.table(Whole_database_GWAS_Phar_Clin,file = "Whole_database_GWAS_Phar_Clin.txt",row.names = FALSE,quote = FALSE,sep = "\t")

## reread
library(data.table)
Whole_database_GWAS_Phar_Clin <- fread("/Users/xiyas/V2_Genome_reporting/Whole_database_GWAS_Phar_Clin.txt",sep ="\t" )

####separate again for mapping sza id for every table
Whole_database_GWAS_Phar_Clin$Type <-sapply(strsplit(Whole_database_GWAS_Phar_Clin$INFO,split = ";"),"[[",1)
Whole_database_GWAS_Phar_Clin$SZAID <-sapply(strsplit(Whole_database_GWAS_Phar_Clin$INFO,split = ";"),"[[",2)

Whole_database_GWAS_Phar_Clin$Type <- gsub(pattern = "Type=",replacement = "",x = Whole_database_GWAS_Phar_Clin$Type)
Whole_database_GWAS_Phar_Clin$SZAID <- gsub(pattern = "SZAID=",replacement = "",x = Whole_database_GWAS_Phar_Clin$SZAID)
write.table(Whole_database_GWAS_Phar_Clin,file = "Whole_database_GWAS_Phar_Clin_for_check.txt",row.names = FALSE,quote = FALSE,sep = "\t")

##afterwards this file is transformed into:Whole_out_sorted.vcf
#cat Whole_database_GWAS_Phar_Clin.txt | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Whole_out_sorted.vcf
#######need to convert it to vcf file before annnotation 
##fileformat=VCFv4.1
##fileDate=2022-08-16
##source=ClinVar
##reference=GRCh38
##ID=<Description="ClinVar Variation ID">
##INFO=<ID=SZAID,Number=.,Type=String,Description="SZAvarID">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

SZAID<- Whole_database_GWAS_Phar_Clin[Whole_database_GWAS_Phar_Clin$Type == "ClinP_LP_var",]
SZAID <- SZAID[,c("ID","SZAID")]
############################

#this generating files for cheng

SZAID<- Whole_database_GWAS_Phar_Clin[Whole_database_GWAS_Phar_Clin$Type == "ClinP_LP_var",]
SZAID <- SZAID[,c("ID","SZAID")]
clin_vcf <- fread("clinvar_20221129.P.LP.modify.vcf.gz")
clin_vcf$ID <- as.character(clin_vcf$ID)
Merged_clinvar_vcf <- merge(clin_vcf,SZAID,by.x ="ID",by.y = "ID")
colnames(Merged_clinvar_vcf)
Merged_clinvar_vcf <- Merged_clinvar_vcf[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","SZAID")]
Merged_clinvar_vcf <- Merged_clinvar_vcf[order(Merged_clinvar_vcf$SZAID),]
#Merged_clinvar_vcf <- Merged_pharma_vcf[,c("CHROM","POS","ID","REF","ALT","Clinical_annotation_ID","Genotype/Allele",
#                                          "Annotation Text","Gene","Level of Evidence","Score","Phenotype Category","Drug(s)","Phenotype(s)","pharma_genotype","SZAID")]
write.table(Merged_clinvar_vcf,file = "Merged_clinvar_vcf_szaid_1204.txt",quote = FALSE,sep = "\t",row.names = FALSE)
