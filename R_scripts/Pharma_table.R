library(data.table)
library(dplyr)

# Old version until 2024.6 =========================
#setwd("/Users/xiyas/V2_Genome_reporting/V2_database_prepare")

setwd('/Users/xiyas/ATPM_Project/database-v3')

# New version start from 2024.8 =========================
clin_ann_alleles <-fread("clinicalAnnotations_20240812/clinical_ann_alleles.tsv",header = TRUE,fill = TRUE,sep = "\t")
clin_annotations <- fread("clinicalAnnotations_20240812/clinical_annotations.tsv",header = TRUE,fill = TRUE,sep = "\t")
clin_ann_evidence <- fread("clinicalAnnotations_20240812/clinical_ann_evidence.tsv",header = TRUE,fill = TRUE,sep = "\t")
colnames(clin_ann_alleles)[1] <- "Clinical_annotation_ID"
colnames(clin_annotations)[1] <- "Clinical_annotation_ID"
Pharma_table <- merge(x=clin_ann_alleles,y=clin_annotations,by = "Clinical_annotation_ID", all.x = TRUE)

##############2024.4 A Haplotype_db is needed to construct: ============================
#This is basically as same as the Merged_Pharma_vcf,
##but no need to merge, just need the original PharmGKB file + complelement the thing previously removed + find available 
# rsID in the /Users/xiyas/V2_Genome_reporting/PharmVar_v6.1.numbers
haplotype_rs_ref <- read.csv("/Users/xiyas/ATPM_Project/database-v3/Alldata0510.csv",header = TRUE,fill = TRUE)
# new method keep as much as possible, both rsSNP and haplotype with rsSNP id
Pharma_table$connected_star <- paste0(Pharma_table$Gene, Pharma_table$`Genotype/Allele`)
# adding rs id
Pharma_table$matched_rs_id <- haplotype_rs_ref$rsID[match(Pharma_table$connected_star,haplotype_rs_ref$Allele)]
# old method only keep the variants with rs ID ======================
Pharma_table_without_hap <- Pharma_table[grepl("rs", Pharma_table$`Variant/Haplotypes`),]
haplotype <- Pharma_table[!grepl("rs", Pharma_table$`Variant/Haplotypes`),]
#get these haplotype list with the format GENE*1, GENE*2, etc
#The paste0() function concatenates strings without any separator.
haplotype$connected_star <- paste0(haplotype$Gene, haplotype$`Genotype/Allele`)
# get a list of how many these haplytype contains rs ID  correspondingly
#haplotype$connected_star %in% haplotype_rs_ref$Allele
# adding rs id
Pharma_table$matched_rs_id <- haplotype_rs_ref$rsID[match(Pharma_table$connected_star,haplotype_rs_ref$Allele)]
Pharma_table$matched_REF <- haplotype_rs_ref$Reference.Allele[match(Pharma_table$connected_star,haplotype_rs_ref$Allele)]
Pharma_table$matched_ALT <- haplotype_rs_ref$Variant.Allele[match(Pharma_table$connected_star,haplotype_rs_ref$Allele)]

haplotype$matched_rs_id <- haplotype_rs_ref$rsID[match(haplotype$connected_star,haplotype_rs_ref$Allele)]
haplotype$matched_REF <- haplotype_rs_ref$Reference.Allele[match(haplotype$connected_star,haplotype_rs_ref$Allele)]
haplotype$matched_ALT <- haplotype_rs_ref$Variant.Allele[match(haplotype$connected_star,haplotype_rs_ref$Allele)]

#table(haplotype$connected_star %in% haplotype_rs_ref$Allele)
#FALSE  TRUE 
#947  1239 

# Summary Variant count by Level of evidence ==========
# it suggests that there are multiple unique Variant/Haplotypes within each group defined by Level of Evidence
Level_of_evidence <- Pharma_table_without_hap %>%
  group_by(`Level of Evidence`) %>%
  summarise(Count = n_distinct(`Variant/Haplotypes`))

# varints that Tier 1 VIP and also strong level of evidence
Tier_1_VIP  %>% group_by(`Level of Evidence`) %>%
  summarise(Count = n_distinct(`Variant/Haplotypes`))

# if just using Tier 1 VIP genes, you get too less variants indentified in 
Target_variants <- Pharma_table_without_hap %>%
  filter(`Level of Evidence` %in% c('1A', '1B', '2A','2B'))
Target_variants_rs_id <- unique(Target_variants$`Variant/Haplotypes`)


# write tables for VEP annotation to get chromosome location ==========
## 2833 normal variants
write.table(unique(Pharma_table_without_hap$`Variant/Haplotypes`),file = "Pharma_table_rs_id_2024.txt",row.names = FALSE,quote = FALSE,na="")
## adding 176 rsID
write.table(unique(haplotype$matched_rs_id),file = "Pharma_table_rs_id_2024.txt",row.names = FALSE,col.names = FALSE, quote = FALSE,append = TRUE,na="")

#sum = 2833+176 = 3009 Pharmaco variants
list = union(unique(Pharma_table_without_hap$`Variant/Haplotypes`),unique(haplotype$matched_rs_id))
# put this Pharma_table_rs_id_2024.txt to VEP Online !!!

#write.table(Pharma_table_without_hap,file = "Pharma_table_without_hap.txt",row.names = FALSE,quote = FALSE)
#colnames(Pharma_table_without_hap)[5] <- "ID"

#using vep website online to get the Pharama_rs_id_vep.vcf ================================
#Pharma_vcf <- fread("Pharama_rs_id_vep.vcf")
Pharma_vcf <- fread("/Users/xiyas/ATPM_Project/database-v3/Pharma_rs_id_vep_2024.vcf")

#Pharma_vcf <- Pharma_vcf[Pharma_vcf$ID %in% Pharma_table_without_hap$ID]
table(Pharma_vcf$ID %in% Pharma_table$`Variant/Haplotypes` | Pharma_vcf$ID %in% Pharma_table$matched_rs_id)
# make the Sum_ID equal to Variant/Haplotypes while it start with rs and equal to matched_rs_id while it not start with rs
Pharma_table$Sum_ID <- ifelse(grepl("^rs", Pharma_table$`Variant/Haplotypes`), 
                              Pharma_table$`Variant/Haplotypes`, 
                             Pharma_table$matched_rs_id)
table(Pharma_vcf$ID %in% Pharma_table$Sum_ID)

### ALL TRUE 

###chcek duplicated row
dup <- Pharma_vcf[duplicated(Pharma_vcf$ID),]
colnames(Pharma_vcf)[1] <- "CHROM"

####is chromosome problems
Pharma_vcf <- Pharma_vcf[Pharma_vcf$CHROM %in% c("MT","X","Y",1:22),]
#check which one in pharmacogenomics table is not annotated by vep
#setdiff(Pharma_table_without_hap$ID,Pharma_vcf$ID)
setdiff(Pharma_table$Sum_ID,Pharma_vcf$ID)
setdiff(Pharma_vcf$ID,Pharma_table$Sum_ID)
##remove the one vep failed to find
#Pharma_table_without_hap<- Pharma_table_without_hap[!which(Pharma_table_without_hap$ID %in% c("rs1799735","rs4630","rs36056065"))]
#
Merged_pharma_vcf <- merge(Pharma_table,Pharma_vcf,by.x="Sum_ID",by.y= "ID",all.y = TRUE)

##adding pharma_genoytpe ====================
###
# Loop through each row in the Merged_pharma_vcf data frame
for (row in 1:nrow(Merged_pharma_vcf)) {
  # Check if the current row's Variant/Haplotypes starts with "rs"
  if (grepl("^rs", Merged_pharma_vcf$`Variant/Haplotypes`[row])) {
    # Construct the pharma_genotype using REF and Genotype/Allele values
    Merged_pharma_vcf$pharma_genotype[row] <- 
      paste0(Merged_pharma_vcf$REF[row], "/", Merged_pharma_vcf$REF[row], ">", 
             paste0(strsplit(Merged_pharma_vcf$`Genotype/Allele`[row], split = "")[[1]], collapse = "/"))
  } else {
    # Construct the pharma_genotype using matched_REF and matched_ALT values
    Merged_pharma_vcf$pharma_genotype[row] <- 
      paste0(Merged_pharma_vcf$matched_REF[row], "/", Merged_pharma_vcf$matched_REF[row], ">", 
             Merged_pharma_vcf$matched_ALT[row], "/", Merged_pharma_vcf$matched_ALT[row])
  }
}
###ADDING GENOTYPES ============================
#Merged_pharma_vcf$risk_genotype <- paste0(Merged_GWAS_vcf$REF,"/",Merged_GWAS_vcf$REF, ">", Merged_GWAS_vcf$REF,"/", Merged_GWAS_vcf$risk_allele)
#Merged_GWAS_vcf$SZAID.x <- NULL
#Merged_GWAS_vcf$SZAID.y <- NULL
########map with sza id
#Whole_database_GWAS_Phar_Clin <- read.delim("~/Report_genome/Whole_database_GWAS_Phar_Clin_for_check.txt")
#SZAID<- Whole_database_GWAS_Phar_Clin[Whole_database_GWAS_Phar_Clin$Type == "Pharma_var",]
#SZAID <- SZAID[,c("ID","SZAID")]
#Merged_pharma_vcf <- merge(Merged_pharma_vcf,SZAID,by.x ="ID",by.y = "ID")


#Merged_pharma_vcf_old <-read.delim('/Users/xiyas/V2_Genome_reporting/database-file-v2/Merged_Pharma_vcf.txt')
write.table(Merged_pharma_vcf,file = "/Users/xiyas/ATPM_Project/database-v3/Merged_Pharma_vcf_2024.txt",quote = FALSE,sep = "\t",row.names = FALSE)


###then this is for generating sza id -----------------------------------
Merged_pharma_vcf_for_sza_id <- Merged_pharma_vcf[,c("CHROM","POS","Sum_ID","REF","ALT","QUAL","FILTER","INFO")]
Merged_pharma_vcf_for_sza_id$INFO <- NULL
###remove duplicateds
Merged_pharma_vcf_for_sza_id <- Merged_pharma_vcf_for_sza_id %>% distinct()
Merged_pharma_vcf_for_sza_id$INFO <- rep("Type=Pharma_var;")
Merged_pharma_vcf_for_sza_id$INFO <- paste0(Merged_pharma_vcf_for_sza_id$INFO,"SZAID=SZAvar")
colnames(Merged_pharma_vcf_for_sza_id)[3] <-"ID" 

### [[[[Summary: 2964 PGx variant]]]]

###mapped with sza id ? ==================
SZAID<- Whole_database_GWAS_Phar_Clin[Whole_database_GWAS_Phar_Clin$Type == "Pharma_var",]
SZAID <- SZAID[,c("ID","SZAID")]
Merged_pharma_vcf <- merge(Merged_pharma_vcf,SZAID,by.x ="ID",by.y = "ID")
Merged_pharma_vcf <- Merged_pharma_vcf[,c("CHROM","POS","ID","REF","ALT","Clinical_annotation_ID","Genotype/Allele",
 "Annotation Text","Gene","Level of Evidence","Score","Phenotype Category","Drug(s)","Phenotype(s)","pharma_genotype","SZAID")]
write.table(Merged_pharma_vcf,file = "Pharma_info.txt",row.names = FALSE,quote = FALSE,sep = "\t")

