##quote =" " is important!!!!!
library(dplyr)
library(tidyr)
setwd('/Users/xiyas/V2_Genome_reporting/V2_database_prepare/')
GWAS_ori <- read.table("gwas_catalog_v1.0.2-associations_e107_r2022-08-21.tsv",header = TRUE,fill = TRUE,sep = "\t",quote = "")
#GWAS_ori <- readLines(file("gwas_catalog_v1.0.2-associations_e107_r2022-08-21.tsv", encoding = "utf-8"))
colnames(GWAS_ori)

GWAS <- GWAS_ori[,c("DISEASE.TRAIT","CHR_ID","CHR_POS","REPORTED.GENE.S.",
                    "MAPPED_GENE","STRONGEST.SNP.RISK.ALLELE","SNPS","OR.or.BETA","PUBMEDID","CONTEXT","P.VALUE","STUDY.ACCESSION","MAPPED_GENE","CONTEXT")]
#remove rows can not use in first version: muiti-loci, and loci with ?
#multi-loci
GWAS <- GWAS[-which(lengths(sapply(strsplit(GWAS$STRONGEST.SNP.RISK.ALLELE, ";"), `[`, simplify=FALSE)) >1,arr.ind = TRUE),]
GWAS <- GWAS[-which(lengths(sapply(strsplit(GWAS$STRONGEST.SNP.RISK.ALLELE, ","), `[`, simplify=FALSE)) >1,arr.ind = TRUE),]
GWAS <- GWAS[-which(lengths(sapply(strsplit(GWAS$STRONGEST.SNP.RISK.ALLELE, " +"), `[`, simplify=FALSE)) >1,arr.ind = TRUE),]

##non-snps, like deletion
GWAS <- GWAS[grepl("rs", GWAS$STRONGEST.SNP.RISK.ALLELE),]

GWAS<-  GWAS%>% separate(STRONGEST.SNP.RISK.ALLELE, c('risk_snp', 'risk_allele'),sep= "-")

#remove risk allele is NA

GWAS <- GWAS[!(is.na(GWAS$risk_allele)),]

#remove risk allele is ?
GWAS <- GWAS[!grepl("\\?", GWAS$risk_allele),]

#write.table(unique(GWAS$risk_snp),file = "GWAS_table_rs_id.txt",row.names = FALSE,quote = FALSE)

#########After annotation this time I run it in uppmax vep
GWAS_vcf<- fread("GWAS_table_rs_id_annotated.txt")
####is chromosome problems
colnames(GWAS_vcf)[1] <- "CHROM"
GWAS_vcf <- GWAS_vcf[GWAS_vcf$CHROM %in% c("MT","X","Y",1:22),]

setdiff(GWAS$risk_snp,GWAS_vcf$ID)
###Remove snps didn't annotated

GWAS<- GWAS[which(GWAS$risk_snp %in% GWAS_vcf$ID),]
colnames(GWAS)[6] <- "ID"
Merged_GWAS_vcf <- merge(GWAS,GWAS_vcf,by.x = "ID",by.y = "ID",all.x = TRUE)


###ADDING GENOTYPES
Merged_GWAS_vcf$risk_genotype <- paste0(Merged_GWAS_vcf$REF,"/",Merged_GWAS_vcf$REF, ">", Merged_GWAS_vcf$REF,"/", Merged_GWAS_vcf$risk_allele)
Merged_GWAS_vcf$SZAID.x <- NULL
Merged_GWAS_vcf$SZAID.y <- NULL
########map with sza id
SZAID<- Whole_database_GWAS_Phar_Clin[Whole_database_GWAS_Phar_Clin$Type == "GWAS_var",]
SZAID <- SZAID[,c("ID","SZAID")]
Merged_GWAS_vcf <- merge(Merged_GWAS_vcf,SZAID,by.x ="ID",by.y = "ID")

###save the file
write.table(Merged_GWAS_vcf,file = "/Users/xiyas/V2_Genome_reporting/Merged_GWAS_vcf.txt",quote = FALSE,sep = "\t",row.names = FALSE)

Merged_GWAS_vcf <- read.delim("/Users/xiyas/V2_Genome_reporting/Merged_GWAS_vcf.txt")
###This is for generating szaids
Merged_GWAS_vcf_for_sza_id <- Merged_GWAS_vcf[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")]
Merged_GWAS_vcf_for_sza_id$INFO <- NULL
###remove duplicateds
Merged_GWAS_vcf_for_sza_id <- Merged_GWAS_vcf_for_sza_id %>% distinct()

Merged_GWAS_vcf_for_sza_id$INFO <- rep("Type=GWAS_var;")
Merged_GWAS_vcf_for_sza_id$INFO <- paste0(Merged_GWAS_vcf_for_sza_id$INFO,"SZAID=SZAvar")


write.table(clin_vcf_sza_id,file = "clin_vcf_sza_id.vcf",quote = FALSE,sep = "\t",row.names = FALSE)



