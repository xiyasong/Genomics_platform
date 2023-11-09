library(data.table)
setwd("~/Report_genome/New_v2_database/")
clin_ann_alleles <- fread("clinicalAnnotations/clinical_ann_alleles.tsv",header = TRUE,fill = TRUE,sep = "\t")
clin_annotations <- fread("clinicalAnnotations/clinical_annotations.tsv",header = TRUE,fill = TRUE,sep = "\t")
clin_ann_evidence <- fread("clinicalAnnotations/clinical_ann_evidence.tsv",header = TRUE,fill = TRUE,sep = "\t")

colnames(clin_ann_alleles)[1] <- "Clinical_annotation_ID"
colnames(clin_annotations)[1] <- "Clinical_annotation_ID"

Pharma_table <- merge(x=clin_ann_alleles,y=clin_annotations,by = "Clinical_annotation_ID", all.x = TRUE)
Pharma_table_without_hap <- Pharma_table[grepl("rs", Pharma_table$`Variant/Haplotypes`),]
write.table(unique(Pharma_table_without_hap$`Variant/Haplotypes`),file = "Pharma_table_rs_id.txt",row.names = FALSE,quote = FALSE)
write.table(Pharma_table_without_hap,file = "Pharma_table_without_hap.txt",row.names = FALSE,quote = FALSE)

colnames(Pharma_table_without_hap)[5] <- "ID"

#using vep website online
Pharma_vcf <- fread("Pharama_rs_id_vep.vcf")
Pharma_vcf <- Pharma_vcf[Pharma_vcf$ID %in% Pharma_table_without_hap$ID]

###chcek duplicated row
Pharma_vcf[duplicated(Pharma_vcf$ID),]
colnames(Pharma_vcf)[1] <- "CHROM"

####is chromosome problems
Pharma_vcf <- Pharma_vcf[Pharma_vcf$CHROM %in% c("MT","X","Y",1:22),]



#check which one in pharmacogenomics table is not annotated by vep
setdiff(Pharma_table_without_hap$ID,Pharma_vcf$ID)

##remove the one vep failed to find
Pharma_table_without_hap<- Pharma_table_without_hap[!which(Pharma_table_without_hap$ID %in% c("rs1799735","rs4630","rs36056065"))]
Merged_pharma_vcf <- merge(Pharma_table_without_hap,Pharma_vcf,by.x = "ID",by.y = "ID",all.x = TRUE)

##adding pharma_genoytpe ???



###
#for (row in 1:nrow(Pharma_vcf)) {
#  Merged_pharma_vcf$pharma_genotype[row] <- 
#    paste0(Merged_pharma_vcf$REF[row],"/",Merged_pharma_vcf$REF[row], ">", 
#           paste0(strsplit(Merged_pharma_vcf$Genotype.Allele[row],split = "")[[1]],collapse = "/"))
#}


###ADDING GENOTYPES
#Merged_pharma_vcf$risk_genotype <- paste0(Merged_GWAS_vcf$REF,"/",Merged_GWAS_vcf$REF, ">", Merged_GWAS_vcf$REF,"/", Merged_GWAS_vcf$risk_allele)
#Merged_GWAS_vcf$SZAID.x <- NULL
#Merged_GWAS_vcf$SZAID.y <- NULL
########map with sza id
#Whole_database_GWAS_Phar_Clin <- read.delim("~/Report_genome/Whole_database_GWAS_Phar_Clin_for_check.txt")
#SZAID<- Whole_database_GWAS_Phar_Clin[Whole_database_GWAS_Phar_Clin$Type == "Pharma_var",]
#SZAID <- SZAID[,c("ID","SZAID")]
#Merged_pharma_vcf <- merge(Merged_pharma_vcf,SZAID,by.x ="ID",by.y = "ID")
Merged_pharma_vcf <-read.delim('/Users/xiyas/V2_Genome_reporting/Merged_Pharma_vcf.txt')
for (row in 1:nrow(Merged_pharma_vcf)) {
  Merged_pharma_vcf$pharma_genotype[row] <- 
    paste0(Merged_pharma_vcf$REF[row],"/",Merged_pharma_vcf$REF[row], ">", 
           paste0(strsplit(Merged_pharma_vcf$Genotype.Allele[row],split = "")[[1]],collapse = "/"))
}

write.table(Merged_pharma_vcf,file = "/Users/xiyas/V2_Genome_reporting/Merged_Pharma_vcf.txt",quote = FALSE,sep = "\t",row.names = FALSE)

###then this is for generating sza id
Merged_pharma_vcf_for_sza_id <- Merged_pharma_vcf[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")]
Merged_pharma_vcf_for_sza_id$INFO <- NULL
###remove duplicateds
Merged_pharma_vcf_for_sza_id <- Merged_pharma_vcf_for_sza_id %>% distinct()
Merged_pharma_vcf_for_sza_id$INFO <- rep("Type=Pharma_var;")
Merged_pharma_vcf_for_sza_id$INFO <- paste0(Merged_pharma_vcf_for_sza_id$INFO,"SZAID=SZAvar")


###mapped with sza id
SZAID<- Whole_database_GWAS_Phar_Clin[Whole_database_GWAS_Phar_Clin$Type == "Pharma_var",]
SZAID <- SZAID[,c("ID","SZAID")]
Merged_pharma_vcf <- merge(Merged_pharma_vcf,SZAID,by.x ="ID",by.y = "ID")
Merged_pharma_vcf <- Merged_pharma_vcf[,c("CHROM","POS","ID","REF","ALT","Clinical_annotation_ID","Genotype/Allele",
 "Annotation Text","Gene","Level of Evidence","Score","Phenotype Category","Drug(s)","Phenotype(s)","pharma_genotype","SZAID")]
write.table(Merged_pharma_vcf,file = "Pharma_info.txt",row.names = FALSE,quote = FALSE,sep = "\t")

