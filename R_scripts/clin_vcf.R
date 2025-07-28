###Results whole analysis 0525-2023

## Gene-centric database checking 
table(GeneDB$Gene.Disease.confidence.level)

## matching case
# Assuming you have a data frame called 'data' with 'genes' and 'diseases' columns

# Count the number of diseases each gene matches
gene_counts <- table(GeneDB$Genes)

# Count the occurrences of genes matching a specific number of diseases
matching_counts <- table(gene_counts)

total_genes <- length(gene_counts)


# Print the results
for (i in seq_along(matching_counts)) {
  num_matches <- names(matching_counts)[i]
  count <- matching_counts[i]
  proportion <- count / total_genes
  
  cat("Number of genes matching", num_matches, "disease(s):", count, "\n")
  cat("Proportion:", proportion, "\n\n")
}

###check the disease inheritances


#########The variants features
library(data.table)
Whole_database_GWAS_Phar_Clin <- fread("/Users/xiyas/V2_Genome_reporting/Whole_database_GWAS_Phar_Clin.txt",sep ="\t" )
Whole_database_GWAS_Phar_Clin$Type <-sapply(strsplit(Whole_database_GWAS_Phar_Clin$INFO,split = ";"),"[[",1)
Whole_database_GWAS_Phar_Clin$SZAID <-sapply(strsplit(Whole_database_GWAS_Phar_Clin$INFO,split = ";"),"[[",2)
table(Whole_database_GWAS_Phar_Clin$Type)


#########clinvar database features
library(tidyverse)
library(data.table)
clin_vcf <- fread("/Users/xiyas/V2_Genome_reporting/clinvar_20221129.P.LP.modify.vcf.gz")

INFO_field <- strsplit(as.character(clin_vcf$INFO),';')
####NEED exact match 
##fixed function is not work, 
CLNSIG<- lapply(INFO_field, grep, pattern = 'CLNSIG=', value = TRUE)
clin_vcf$CLNSIG =  CLNSIG
clin_vcf$CLNSIG = gsub("CLNSIG=","",clin_vcf$CLNSIG)
clin_vcf$CLNSIG <- unlist(clin_vcf$CLNSIG)
clin_vcf$CLNSIG <- unlist(lapply(clin_vcf$CLNSIG, function(x) strsplit(x, "\\|")[[1]][1]))
clin_vcf$CLNSIG <- unlist(lapply(clin_vcf$CLNSIG, function(x) strsplit(x, "\\/")[[1]][1]))
clin_vcf$CLNSIG = gsub("Likely_pathogenic_low_penetrance","Likely_pathogenic",clin_vcf$CLNSIG)
table(clin_vcf$CLNSIG)
head(clin_vcf$CLNSIG)

Gene_Info<- lapply(INFO_field, grep, pattern = 'GENEINFO=', value = TRUE)
clin_vcf$Gene_Info =  Gene_Info
clin_vcf$Gene_Info = gsub("GENEINFO=","",clin_vcf$Gene_Info)
clin_vcf$Gene_Info <- unlist(lapply(clin_vcf$Gene_Info, function(x) strsplit(x, ":")[[1]][1]))

Rev_stat<- lapply(INFO_field, grep, pattern = 'REVSTAT=', value = TRUE)
clin_vcf$Rev_stat =  Rev_stat
clin_vcf$Rev_stat = gsub("CLNREVSTAT=","",clin_vcf$Rev_stat)
#clin_vcf$Gene_Info <- unlist(lapply(clin_vcf$Gene_Info, function(x) strsplit(x, ":")[[1]][1]))

var_counts <- count(clin_vcf,Gene_Info)

top_var <- var_counts %>%
  arrange(desc(n)) %>%
  head(50)


counted_df <- clin_vcf %>% group_by(Rev_stat, CLNSIG)%>% summarise(count = n())
counted_df <- counted_df %>%
  group_by(CLNSIG,Rev_stat) %>%
  mutate(proportion = count / sum(count))
custom_labels <- c("Conflicts", "Likely pathogenic", "Pathogenic")
# Create a grouped bar plot with CLNSIG on the x-axis and fill colors representing Rev_stat
ggplot(counted_df, aes(x = CLNSIG, y = count, fill = Rev_stat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Count of Review status within groups of Clinvar pathogenicity", x = "CLNSIG", y = "Counts")+
  xlab("Pathogenicity") +theme_bw() +
  scale_fill_jco() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid = element_blank())+ 
  scale_x_discrete(labels = custom_labels)
# Rotate x-axis labels for better readability


#####Pharmacogenetics variants

#pharma_C <- table(Merged_pharma_vcf$Drug.s.)
#pharma_C <- as.data.frame(pharma_C)
#colnames(pharma_C) <- c("Value", "Frequency")
#pharma_C <- pharma_C[order(pharma_C$Frequency, decreasing = TRUE), ]

#The unique records (unique variants)

unique_rows <- Merged_pharma_vcf[!duplicated(Merged_pharma_vcf$ID),]
pharma_C <- table(unique_rows$Drug.s.)
pharma_C <- as.data.frame(pharma_C)
colnames(pharma_C) <- c("Value", "Frequency")
pharma_C <- pharma_C[order(pharma_C$Frequency, decreasing = TRUE), ]

##### GWAS variants 
#GWAS_C <- table(Merged_GWAS_vcf$DISEASE.TRAIT)
#pharma_C <- as.data.frame(pharma_C)
#colnames(pharma_C) <- c("Value", "Frequency")
#pharma_C <- pharma_C[order(pharma_C$Frequency, decreasing = TRUE), ]
#The unique records (unique variants)

unique_rows <- Merged_GWAS_vcf[!duplicated(Merged_GWAS_vcf$ID),]

GWAS_C <- table(unique_rows$DISEASE.TRAIT)
GWAS_C <- as.data.frame(GWAS_C)
colnames(GWAS_C) <- c("Value", "Frequency")
GWAS_C <- GWAS_C[order(GWAS_C$Frequency, decreasing = TRUE), ]

