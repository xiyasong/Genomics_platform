####load the package-----------
library(dplyr)
library(ggsci)
library(ggpubr)
library(ggtext)
library(tidyverse)
library(ggsci)
library(reshape2)
library(ggplot2)
library(UpSetR)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("krassowski/complex-upset")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
# Get the summary DB combines GeneDB and DiseaseDB and Inheritance ----------
DiseaseDB <- read.delim2("/Users/xiyas/V2_Genome_reporting/database-file-v2/DiseaseDB_version_2_Beta_curated.txt")
GeneDB <- read.delim2("/Users/xiyas/V2_Genome_reporting/database-file-v2/GeneDB.txt")
GeneDB_new <- read.delim2("/Users/xiyas/ATPM_Project/database-v3/GeneDB_GenCC.txt")
OMIM_inheritance_DB <- read.delim2('/Users/xiyas/V2_Genome_reporting/database-file-v2/pheno_OMIM_all.txt') 
OMIM_inheritance_DB$phenotypeMimNumber <- as.character(OMIM_inheritance_DB$phenotypeMimNumber)
###Mapping database
###left join, only kept the diseases that in the GeneDB  
Summary_DB <- left_join(GeneDB,DiseaseDB,by= "SZAdiseaseID")
Summary_DB <- left_join(Summary_DB,OMIM_inheritance_DB,by= c("DiseaseMIM"="phenotypeMimNumber"))
Summary_DB <- Summary_DB %>% filter(!duplicated(cbind(Genes,SZAdiseaseID,Gene.Disease.confidence.level,Target.group)))

AD_genes <- 
  unique(
    Summary_DB %>% filter(grepl("Autosomal dominant",inheritances))%>% 
      .$Genes
  )
AR_genes <- 
  unique(
    Summary_DB %>% filter(grepl("Autosomal recessive",inheritances))%>% 
      .$Genes
  )

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
# Figure S2A ==============
# Define custom colors
custom_colors  <- pal_lancet("lanonc", alpha = 0.9)(6)
# Group by genes and count the unique diseases associated with each gene in Summary_DB =======
gene_counts <- Summary_DB %>%
  group_by(Genes) %>%
  summarise(NumDiseases = n_distinct(SZAdiseaseID))

# Now, you can calculate the number of genes with only one disease and two or more diseases =====
num_genes_one_disease <- sum(gene_counts$NumDiseases == 1)
num_genes_two_diseases <- sum(gene_counts$NumDiseases == 2)
num_genes_three_diseases <- sum(gene_counts$NumDiseases == 3)
num_genes_more_diseases <- sum(gene_counts$NumDiseases > 3)
counts <- c(
  "Matching one Disease" = num_genes_one_disease,
  "Matching two Diseases" = num_genes_two_diseases,
  "Matching three Diseases" = num_genes_three_diseases,
  "Matching more than three Diseases" = num_genes_more_diseases
)
# Calculate proportions
proportions <- round(counts / sum(counts) * 100, 1)
# Define custom colors
#custom_colors <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73")
# Create a pie chart with proportions labeled and custom colors
pdf("~/V2_Genome_reporting/R_scripts/figures/Figure_S2A.pdf", width = 8.5, height = 4)
pie(counts, labels = paste(names(counts), "(", proportions, "%)"), 
    main = "Distribution of Genes by Number of Diseases", 
    col = custom_colors,cex = 1)
dev.off()
# ----------------------
disease_couns <- Summary_DB %>%
  group_by(Disease) %>%
  summarise(NumGenes = n_distinct(Genes))

num_genes_one_disease <- sum(disease_couns$NumGenes == 1)
num_genes_two_diseases <- sum(disease_couns$NumGenes == 2)
num_genes_three_diseases <- sum(disease_couns$NumGenes == 3)
num_genes_more_diseases <- sum(disease_couns$NumGenes > 3)
counts <- c(
  "Matching one Gene" = num_genes_one_disease,
  "Matching two Genes" = num_genes_two_diseases,
  "Matching Three Genes" = num_genes_three_diseases,
  "Matching more than three Genes" = num_genes_more_diseases
)
# Calculate proportions
proportions <- round(counts / sum(counts) * 100, 1)
pdf("~/V2_Genome_reporting/R_scripts/figures/Figure_S2B.pdf", width = 8.5, height = 4)

pie(counts, labels = paste(names(counts), "(", proportions, "%)"), 
    main = "Distribution of Diseases by Number of Genes", 
    col = custom_colors,cex = 1)
dev.off()

#=Figure S2C ======================
Summary_DB$inheritances <-sapply(strsplit(Summary_DB$inheritances,","), `[`, 1)
# Sample data
table(is.na(Summary_DB$inheritances))
data <- c("Autosomal dominant" =2475, "Autosomal recessive"=3535, 
          "X-linked dominant"=73, "X-linked recessive"=217, 
          "no patterns"=3551,"Others"=297)
# Calculate proportions
proportions <- round(data / sum(data) * 100, 1)

#label_pos <- 3.2
# Create a pie chart with proportions labeled and custom colors
pdf("~/V2_Genome_reporting/R_scripts/figures/Figure_S2C.pdf", width = 8.5, height = 4)
pie(data, labels = paste(names(data), "(", proportions, "%)"), 
    main = "Distribution of Gene-disease pairs by inheritance patterns", 
    col = custom_colors,cex = 1)
dev.off()

### Figure S2D ---------------
x <- list(
  set1 <- unique(GeneDB[GeneDB$Target.group == "Carrier-screening",]$Genes),
  set2 <- unique(GeneDB[GeneDB$Target.group == "Newborn-screening",]$Genes),
  set3 <- unique(GeneDB[GeneDB$Target.group == "Health predipositions/Disease risk",]$Genes),
  set4 <- unique(GeneDB[GeneDB$Target.group == "Heriditary-cancer risk syndrome",]$Genes)
)

display_venn(
  x,
  category.names = c("Carrier-screening" , "Newborn-screening" , "Health predipositions", "Heriditary-cancer risk syndrome"),
  fill = pal_jco("default")(4),cat.cex = 2, # 设置类别标签字体大小
  cex = 1,margin = 0.2,cat.dist = 0.15
  
)
### Figure S2E ---------------
x <- list(
  set1 <- unique(GeneDB[GeneDB$Target.group == "Carrier-screening",]$Genes),
  set2 <- unique(GeneDB[GeneDB$Target.group == "Newborn-screening",]$Genes),
  set3 <- unique(GeneDB[GeneDB$Target.group == "Expanded_clinvar_mono_diseases",]$Genes),
  set4 <- unique(GeneDB[GeneDB$Target.group == "Expanded_mono_rare_diseases",]$Genes)
)

display_venn(
  x,
  category.names = c("Carrier-screening" , "Newborn-screening " , "Clinvar genes", "mono_rare"),
  fill = pal_jco("default")(4),cat.cex = 2, # 设置类别标签字体大小
  cex = 2,margin = 0.2,cat.dist = 0.15
)

x <- list(
  set1 <- unique(GeneDB[GeneDB$Target.group == "Carrier-screening",]$SZAdiseaseID),
  set2 <- unique(GeneDB[GeneDB$Target.group == "Newborn-screening",]$SZAdiseaseID),
  set3 <- unique(GeneDB[GeneDB$Target.group == "Health predipositions/Disease risk",]$SZAdiseaseID),
  set4 <- unique(GeneDB[GeneDB$Target.group == "Heriditary-cancer risk syndrome",]$SZAdiseaseID)
)
display_venn(
  x,
  category.names = c("Carrier-screening" , "Newborn-screening " , "Health predipositions", "Heriditary-cancer risk syndrome"),
  fill =  pal_nejm("default")(4)
)
#########UpSet plot
##Figure S2F ---------
x <- list(
  set1 <- unique(GeneDB[GeneDB$Target.group == "Carrier-screening",]$Genes),
  set2 <- unique(GeneDB[GeneDB$Target.group == "Newborn-screening",]$Genes),
  set3 <- unique(GeneDB[GeneDB$Target.group == "Health predipositions/Disease risk",]$Genes),
  set4 <- unique(GeneDB[GeneDB$Target.group == "Heriditary-cancer risk syndrome",]$Genes),
  set5 <- unique(GeneDB[GeneDB$Target.group == "Expanded_clinvar_mono_diseases",]$Genes),
  set6 <- unique(GeneDB[GeneDB$Target.group == "Expanded_mono_rare_diseases",]$Genes)
)

# 生成 named list
names(x) <- c(
  "Carrier-screening", 
  "Newborn-screening",
  "Health predipositions",
  "Heriditary-cancer risk syndrome",
  "Expanded_clinvar_mono_diseases",
  "Expanded_mono_rare_diseases"
)

# 转换为 UpSetR 格式
data_v2 <- UpSetR::fromList(x)

# 绘制 UpSet plot
UpSetR::upset(data_v2, 
              sets = names(x), sets.bar.color = "steelblue",  # 设置集合条的颜色
              matrix.color = "indianred3",       # 设置矩阵点的颜色
              main.bar.color = "darkgreen",  # 设置主条的颜色
              order.by = "freq",             # 按频率排序
              number.angles = 30,            # 数字倾斜角度
              text.scale = 1.5              )
##Figure S2G ---------
x <- list(
  set1 <- unique(GeneDB[GeneDB$Target.group == "Carrier-screening",]$SZAdiseaseID),
  set2 <- unique(GeneDB[GeneDB$Target.group == "Newborn-screening",]$SZAdiseaseID),
  set3 <- unique(GeneDB[GeneDB$Target.group == "Health predipositions/Disease risk",]$SZAdiseaseID),
  set4 <- unique(GeneDB[GeneDB$Target.group == "Heriditary-cancer risk syndrome",]$SZAdiseaseID),
  set5 <- unique(GeneDB[GeneDB$Target.group == "Expanded_clinvar_mono_diseases",]$SZAdiseaseID),
  set6 <- unique(GeneDB[GeneDB$Target.group == "Expanded_mono_rare_diseases",]$SZAdiseaseID)
)
# 生成 named list
names(x) <- c(
  "Carrier-screening", 
  "Newborn-screening",
  "Health predipositions",
  "Heriditary-cancer risk syndrome",
  "Expanded_clinvar_mono_diseases",
  "Expanded_mono_rare_diseases"
)

# 转换为 UpSetR 格式
data_v2 <- UpSetR::fromList(x)
# 绘制 UpSet plot
UpSetR::upset(data_v2, 
              sets = names(x), sets.bar.color = "steelblue",  # 设置集合条的颜色
              matrix.color = "indianred3",       # 设置矩阵点的颜色
              main.bar.color = "darkgreen",  # 设置主条的颜色
              order.by = "freq",             # 按频率排序
              number.angles = 30,            # 数字倾斜角度
              text.scale = 1.5              )


###still pie chart, showing the clinvar variants property
table(clin_vcf$CLNSIG)
data <- c("Pathogenic" =2475, "Autosomal recessive"=3535, 
          "X-linked dominant"=73, "X-linked recessive"=217, 
          "no patterns"=3551,"Others"=297)


Gene_count <- GeneDB%>%group_by(Target.group) %>%
  summarise(NumGenes = n_distinct(Genes))

Disease_count <- GeneDB%>%group_by(Target.group) %>%
  summarise(NumDiseases = n_distinct(SZAdiseaseID))

Total_gene <-length(unique(GeneDB$Genes))
Total_disease <-length(unique(GeneDB$Disease))
Expanded_carrier <-GeneDB %>%filter(Target.group =='Carrier-screening' | Target.group =='Newborn-screening')
Total_carrier_gene <-length(unique(Expanded_carrier$Genes))
Total_carrier_disease <-length(unique(Expanded_carrier$SZAdiseaseID))

#########The variants db features -----------------------
library(data.table)
Whole_database_GWAS_Phar_Clin <- fread("/Users/xiyas/ATPM_Project/database-v3/Whole_database_GWAS_Phar_Clin_2024.txt",sep ="\t" )
Whole_database_GWAS_Phar_Clin$Type <-sapply(strsplit(Whole_database_GWAS_Phar_Clin$INFO,split = ";"),"[[",1)
Whole_database_GWAS_Phar_Clin$SZAID <-sapply(strsplit(Whole_database_GWAS_Phar_Clin$INFO,split = ";"),"[[",2)
table(Whole_database_GWAS_Phar_Clin$Type)
dim(Whole_database_GWAS_Phar_Clin)

Whole_PGx_var1 <-fread('/Users/xiyas/ATPM_Project/database-v3/PGx_annotation_all.levels.drugs.hap250220.txt',sep ="\t")
Whole_PGx_var2 <- fread('/Users/xiyas/ATPM_Project/database-v3/Pharma_table_rs_id_2024.txt',sep ="\t",header = FALSE)

rs_id_set_1 <- unique(Whole_PGx_var1$rsID)
rs_id_set_2 <- unique(Whole_PGx_var2$V1)
length(rs_id_set_1)
#[1] 3115
length(rs_id_set_2)
#[1] 2970
length(unique(union(rs_id_set_1,rs_id_set_2)))
#[1] 3133
all_unique_rsIDs <- unique(union(rs_id_set_1, rs_id_set_2))
# Write rsIDs to a txt file
rs_ids_filtered <- grep("^rs", all_unique_rsIDs, value = TRUE)
length(rs_ids_filtered)
#3106

write.table(rs_ids_filtered,
            file = "/Users/xiyas/ATPM_Project/PGx_all_db_features/all_unique_rsIDs.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# New version start from 2024.8 =========================
setwd('/Users/xiyas/ATPM_Project/database-v3')
clin_ann_alleles <-fread("clinicalAnnotations_20240812/clinical_ann_alleles.tsv",header = TRUE,fill = TRUE,sep = "\t")
clin_annotations <- fread("clinicalAnnotations_20240812/clinical_annotations.tsv",header = TRUE,fill = TRUE,sep = "\t")
clin_ann_evidence <- fread("clinicalAnnotations_20240812/clinical_ann_evidence.tsv",header = TRUE,fill = TRUE,sep = "\t")
colnames(clin_ann_alleles)[1] <- "Clinical_annotation_ID"
colnames(clin_annotations)[1] <- "Clinical_annotation_ID"
Pharma_table <- merge(x=clin_ann_alleles,y=clin_annotations,by = "Clinical_annotation_ID", all.x = TRUE)
# Assuming your column is named 'drug' and the data frame is 'Pharma_table'

# Step 1: Split all drug names by semicolon and unlist into one long vector
all_drugs <- unlist(strsplit(Pharma_table$`Drug(s)`, ";"))
# Step 2: Trim whitespace (if any)
all_drugs <- trimws(all_drugs)
# Step 3: Get unique drug names
unique_drugs <- unique(all_drugs)

# Write unique drug names to a txt file
write.table(unique_drugs,
            file = "/Users/xiyas/ATPM_Project/PGx_all_db_features/all_unique_drug_names.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Step 4: Count how many unique drug names
length(unique_drugs)

length(unique(Pharma_table$`Drug(s)`))
#setdiff(Whole_PGx_var1$rsID,Whole_PGx_var2$V1)

# GWAS is precise, Clinvar should further count CVar for filtered CPLP and Conflicting variants, and PGx should rely on new star allele included db
#Type=ClinP_LP_var     Type=GWAS_var   Type=Pharma_var 
#287187            160618              2964 
# Trait db features
Trait_db <- fread('/Users/xiyas/ATPM_Project/database-v3/Reports_genome_databases_traits_merged_2.txt',sep ="\t",header = TRUE)
length(unique(Trait_db$variants))

#########clinvar database features
library(tidyverse)
library(data.table)
library(dplyr)
lines <- readLines("/Users/xiyas/ATPM_Project/database-v3/clinvar_20240611_PLPC_new_CPLP.vcf.gz")
# Filter out lines starting with ##
filtered_lines <- lines[!grepl("^##", lines)]
clin_vcf <- read.delim2(text = paste(filtered_lines, collapse = "\n"),
                               col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),header = FALSE)
clin_vcf <- clin_vcf%>% distinct()
INFO_field <- strsplit(as.character(clin_vcf$INFO),';')
####NEED exact match 
##fixed function is not work, 
CLNSIG<- lapply(INFO_field, grep, pattern = 'CLNSIG=', value = TRUE)
rsID <- lapply(INFO_field, grep, pattern = 'RS=', value = TRUE)
clin_vcf$rsID <- rsID
clin_vcf$CLNSIG =  CLNSIG
clin_vcf$CLNSIG = gsub("CLNSIG=","",clin_vcf$CLNSIG)
clin_vcf$CLNSIG <- unlist(clin_vcf$CLNSIG)
clin_vcf$CLNSIG <- unlist(lapply(clin_vcf$CLNSIG, function(x) strsplit(x, "\\|")[[1]][1]))
clin_vcf$CLNSIG <- unlist(lapply(clin_vcf$CLNSIG, function(x) strsplit(x, "\\/")[[1]][1]))
clin_vcf$CLNSIG = gsub("Likely_pathogenic,_low_penetrance","Likely_pathogenic",clin_vcf$CLNSIG)
CLNSIGCONF <- lapply(INFO_field, function(x) {
  match <- grep("^CLNSIGCONF=", x, value = TRUE)
  if (length(match) > 0) {
    sub("CLNSIGCONF=", "", match)
  } else {
    NA
  }
})
clin_vcf$CLNSIGCONF <- unlist(CLNSIGCONF)
table(clin_vcf$CLNSIG)
head(clin_vcf$CLNSIG)
dim(clin_vcf)
#287188 #[1] 286274   is the final number !!!!!! 
# 条件 1: CLNSIGCONF 包含 "Pathogenic" 或 "Likely_pathogenic"，且不包含 "benign"
pathogenic_mask <- grepl("Pathogenic|Likely_pathogenic", clin_vcf$CLNSIGCONF, ignore.case = TRUE) &
  !grepl("benign", clin_vcf$CLNSIGCONF, ignore.case = TRUE)

# 条件 2: 如果 CLNSIGCONF 为空，则检查 CLNSIG 列，排除包含 "benign" 或 "likely_benign" 的行
clin_sig_mask <- !grepl("benign|likely_benign", clin_vcf$CLNSIG, ignore.case = TRUE)
# 合并条件：
# - 如果 CLNSIGCONF 满足条件 1，保留；
# - 如果 CLNSIGCONF 为空，且 CLNSIG 满足条件 2，保留。
# final using Clinvar variants is this 
filtered_clin_vcf <- clin_vcf[pathogenic_mask | (is.na(clin_vcf$CLNSIGCONF) & clin_sig_mask), ]

# 2. 统计剩余的各类变异数量
dim(filtered_clin_vcf)
table(clin_vcf$CLNSIG)
table(filtered_clin_vcf$CLNSIG)
#[1] 286274   is the final number !!!!!! 


INFO_field <- strsplit(as.character(filtered_clin_vcf$INFO),';')
Gene_Info<- lapply(INFO_field, grep, pattern = 'GENEINFO=', value = TRUE)
filtered_clin_vcf$Gene_Info =  Gene_Info
filtered_clin_vcf$Gene_Info = gsub("GENEINFO=","",filtered_clin_vcf$Gene_Info)
filtered_clin_vcf$Gene_Info <- unlist(lapply(filtered_clin_vcf$Gene_Info, function(x) strsplit(x, ":")[[1]][1]))

#### supplementary Figure 1
Rev_stat<- lapply(INFO_field, grep, pattern = 'REVSTAT=', value = TRUE)
filtered_clin_vcf$Rev_stat =  Rev_stat
filtered_clin_vcf$Rev_stat = gsub("CLNREVSTAT=","",filtered_clin_vcf$Rev_stat)
#clin_vcf$Gene_Info <- unlist(lapply(clin_vcf$Gene_Info, function(x) strsplit(x, ":")[[1]][1]))

# Function to calculate review stars
calculate_review_star <- function(review) {
  # Handle NA or NULL values
  if (is.na(review) || is.null(review)) {
    review <- ""
  }
  
  # Standardize the review status string
  r <- tolower(review)
  r <- gsub(" ", "", r)
  r <- gsub("_", "", r)
  r <- gsub("&", "", r)
  r <- gsub(",", "", r)
  #print(r)
  
  # Apply review star logic
  if (r %in% c("criteriaprovidedconflictinginterpretations", 
               "criteriaprovidedconflictingclassifications", 
               "criteriaprovidedsinglesubmitter")) {
    return(1)
  } else if (r == "criteriaprovidedmultiplesubmittersnoconflicts") {
    return(2)
  } else if (r == "reviewedbyexpertpanel") {
    return(3)
  } else if (r == "practiceguideline") {
    return(4)
  } else {
    return(0)
  }
}

# Apply review star calculation
filtered_clin_vcf$ReviewStar <- sapply(filtered_clin_vcf$Rev_stat, calculate_review_star)

# Group and count by ReviewStar and CLNSIG
counted_df <- filtered_clin_vcf %>% 
  group_by(ReviewStar, CLNSIG) %>% 
  summarise(count = n()) %>%
  group_by(CLNSIG) %>%
  mutate(proportion = count / sum(count))

# Custom labels for CLNSIG
custom_labels <- c("Conflicts(P/LP)", "Likely pathogenic", "Pathogenic")

# Create plot with ReviewStar
plot <- ggplot(counted_df, aes(x = CLNSIG, y = count, fill = factor(ReviewStar))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Count of Review Stars within groups of Clinvar pathogenicity", 
    x = "Pathogenicity", 
    y = "Counts",
    fill = "Review Star"
  ) +
  theme_bw() +
  scale_fill_nejm() +
  theme(panel.grid = element_blank()) + 
  scale_x_discrete(labels = custom_labels)

# Print or save the plot
print(plot)

### supplementary Table 2  ---------------------------------
library(dplyr)
library(tidyr)

# 1. Reshape GeneDB_new into long format for classification sources
GeneDB_long_classifications <- GeneDB_new %>%
  select(Gene.Disease.confidence.level, starts_with("classification_")) %>%
  pivot_longer(
    cols = starts_with("classification_"),
    names_to = "GenCC_Submitter",
    values_to = "GenCC_Classification"
  )

# 2. Convert NA classifications to string "NA" for counting
GeneDB_long_classifications$GenCC_Classification <- ifelse(
  is.na(GeneDB_long_classifications$GenCC_Classification),
  "NA",
  GeneDB_long_classifications$GenCC_Classification
)

# 3. Count classification values within each confidence level and submitter
GeneDB_classification_counts <- GeneDB_long_classifications %>%
  group_by(Gene.Disease.confidence.level, GenCC_Submitter, GenCC_Classification) %>%
  summarise(Classification_Count = n(), .groups = "drop")

# 4. Calculate proportions within each confidence level + submitter group
GeneDB_classification_summary <- GeneDB_classification_counts %>%
  group_by(Gene.Disease.confidence.level, GenCC_Submitter) %>%
  mutate(
    Total_Pairs = sum(Classification_Count),
    Proportion_Percent = round(Classification_Count / Total_Pairs * 100, 1)
  ) %>%
  ungroup()

# 5. Optional: sort the output
  arrange(Gene.Disease.confidence.level, GenCC_Submitter, desc(Classification_Count))

# View the summary
print(GeneDB_classification_summary)

### much more summary plot

# 1. Create a helper: which rows in GeneDB_new have any GenCC classification?
GenCC_columns <- grep("^classification_", colnames(GeneDB_new), value = TRUE)

GeneDB_new <- GeneDB_new %>%
  mutate(
    Any_GenCC_Classification = apply(select(., all_of(GenCC_columns)), 1, function(row) {
      any(!is.na(row))
    })
  )

# 2. Flatten all GenCC classification values across all sources
GeneDB_long_all_classifications <- GeneDB_new %>%
  select(Gene.Disease.confidence.level, all_of(GenCC_columns)) %>%
  pivot_longer(
    cols = all_of(GenCC_columns),
    names_to = "GenCC_Submitter",
    values_to = "GenCC_Classification"
  ) %>%
  mutate(
    GenCC_Classification = ifelse(is.na(GenCC_Classification), "NA", GenCC_Classification)
  )

# 3. Summarize by confidence level and classification type
ConfidenceLevel_summary <- GeneDB_long_all_classifications %>%
  group_by(Gene.Disease.confidence.level, GenCC_Classification) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  group_by(Gene.Disease.confidence.level) %>%
  mutate(
    Total = sum(Count),
    Proportion = round(Count / Total * 100, 1)
  ) %>%
  arrange(Gene.Disease.confidence.level, desc(Count))

# 4. Optional: View number of gene-disease pairs with/without GenCC classification per confidence level
ConfidenceLevel_pair_summary <- GeneDB_new %>%
  group_by(Gene.Disease.confidence.level) %>%
  summarise(
    Total_GeneDisease_Pairs = n(),
    With_GenCC = sum(Any_GenCC_Classification),
    Without_GenCC = sum(!Any_GenCC_Classification),
    Proportion_With_GenCC = round(With_GenCC / Total_GeneDisease_Pairs * 100, 1)
  )

# View the summaries
print("Summary of classifications by confidence level:")
print(ConfidenceLevel_summary)

print("Summary of gene-disease pairs with/without GenCC classifications:")
print(ConfidenceLevel_pair_summary)
library(writexl)
write_xlsx(
  list(
    "Per_Submission_Summary" = classification_summary_by_confidence,
    "Confidence_level_Summary" = ConfidenceLevel_summary,
    "General" = ConfidenceLevel_pair_summary
  ),
  "/Users/xiyas/ATPM_Project/database-v3/Supplementary_Table_2_GenCC_Summary_Full.xlsx"
)



