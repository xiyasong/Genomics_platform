# get the allele frequency of swe
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(openxlsx)


merge_table_swe <- fread(file= '/Users/xiyas/V2_Genome_reporting/data/merge_table_swe.txt',quote = FALSE,sep = '\t')
merge_table_tur <- fread(file= '/Users/xiyas/V2_Genome_reporting/data/merge_table.txt',quote = FALSE,sep = '\t')


#Pharma_table_without_hap$key <-paste(Pharma_table_without_hap$`Variant/Haplotypes`, Tier_1_VIP$`Genotype/Allele`, sep = "_")
# Database included XX PGx varaints
# Tier_1_VIP ========= first, filter our target variants list to only VIP genes
Tier_1_VIP <- Pharma_table_without_hap %>%
  filter(grepl("Tier 1 VIP", `Level Modifiers`))
length(unique(Tier_1_VIP$`Variant/Haplotypes`))

### 455 variant sites (rs ID)
length(unique(Tier_1_VIP$`Variant/Haplotypes`))


## change genotype to alt allele format
# If Genotype/Allele is G, extracted_char will be G.
# If Genotype/Allele is GC, extracted_char will be C
Tier_1_VIP$extracted_alt <- substring(Tier_1_VIP$`Genotype/Allele`, 
                                       nchar(Tier_1_VIP$`Genotype/Allele`), nchar(Tier_1_VIP$`Genotype/Allele`))

Tier_1_VIP$key <- paste(Tier_1_VIP$`Variant/Haplotypes`, Tier_1_VIP$extracted_alt, sep = "_")

Tier_1_VIP$`Phenotype Category` <- sapply(strsplit(Tier_1_VIP$`Phenotype Category`,";"), `[`, 1)

##shouldn't include 4
Target_variants <- Tier_1_VIP %>%
  filter(`Level of Evidence` %in% c('1A', '1B', '2A','3'))
# 444 target variants
length(unique(Target_variants$`Variant/Haplotypes`))
#Target_variants <-Pharma_table_without_hap %>%
#  filter(`Level of Evidence` %in% c('1A', '1B', '2A'))
Target_variants_rs_id <- unique(Target_variants$`Variant/Haplotypes`)

# adding row annotation ============
# For 'Level of Evidence' (assuming it's a factor with levels like 'High', 'Medium', 'Low')
level_colors <- brewer.pal(n = length(levels(annotation_row$level_of_evidence)), name = "YlOrBr")

# For 'Phenotype Category' (assuming it's also a factor; adjust the number of levels and palette as needed) 
phenotype_colors <- brewer.pal(n = length(levels(annotation_row$phenotype_category)), name = "Paired")

gene_colors = c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",
                "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",
                "orange", "tomato", "coral2", "palevioletred", "violetred", "red2",
                "springgreen2", "yellowgreen", "palegreen4",
                "wheat2", "tan", "tan2", "tan3", "brown",
                "grey70")
annotation_colors <- data.frame(genes = setNames(gene_colors, unique(annotation_row$genes)))

#Target_variants <- Pharma_table_without_hap %>%
#  filter(`Level of Evidence` %in% c('1A', '1B', '2A','2B'))

# start building heatmap ============
Target_variants_rs_id <- unique(Target_variants$`Variant/Haplotypes`)
#Target_variants_rs_id <- unique(Tier_1_VIP$`Variant/Haplotypes`)
#writeLines(Target_variants_rs_id, "rsIDs.txt")


##These numbers are number of biallelic sites
filtered_merge_table_swe <- merge_table_swe[rs_id %in% Target_variants_rs_id]
filtered_merge_table_tur <- merge_table_tur[rs_id %in% Target_variants_rs_id]


###########VEP to retrive the AF data from 1000G and Gnomad =================
##PGx-VIP: rsid with different alt alleles information
kGAF <- fread('/Users/xiyas/V2_Genome_reporting/data//PGx-VIP.txt',sep = '\t',header = TRUE)
kGAF <- as.data.frame(kGAF)
#kGAF <- kGAF[!duplicated(kGAF$`#Uploaded_variation`), ]---this is wrong, should keep all alternative alleles
colnames(kGAF)[1] <- 'Uploaded_variation'

# Merge kGAF with filtered_merge_table_tur to add tur_AF
kGAF$key <- paste(kGAF$Uploaded_variation, kGAF$Allele, sep = "_")
filtered_merge_table_tur$key <- paste(filtered_merge_table_tur$rs_id, filtered_merge_table_tur$ALT, sep = "_")
kGAF <- merge(kGAF, filtered_merge_table_tur[, c("key", "ALT_AF")], by = "key", all.x = TRUE)
colnames(kGAF)[colnames(kGAF) == "ALT_AF"] <- "ALT_AF_TR"

# Merge kGAF with filtered_merge_table_swe to add swe_AF
filtered_merge_table_swe$key <- paste(filtered_merge_table_swe$rs_id, filtered_merge_table_swe$ALT, sep = "_")
kGAF <- merge(kGAF, filtered_merge_table_swe[, c("key", "ALT_AF")], by = "key", all.x = TRUE)
colnames(kGAF)[colnames(kGAF) == "ALT_AF"] <- "ALT_AF_SW"

# I want to add the info from PharmaGKB

#rs9332377_T
kGAF_with_drug <- kGAF %>% left_join(Target_variants, by = "key")


# Remove variants where TR and SW all NA value for allele frequency
# 91 variants allele left
#TR for example, rs72728438 is rs72728438_C;

# several rare PGx variants were only identified in TR, not in SW 1000G and Gnomad v3(Supplementary Table 2).here is as kGAF_only_TR
kGAF <- kGAF[!duplicated(kGAF$key), ] # now start to remove duplicates

annotation_row <- data.frame(
  genes = sapply(strsplit(kGAF$SYMBOL, "-"), `[`, 1)
)
row.names(annotation_row) <- kGAF$key 

kGAF <- kGAF %>% filter(!(is.na(ALT_AF_TR)))
kGAF_only_TR <- kGAF

# Selecting columns by name that contain "AF"
columns_with_AF <- grep("AF", names(kGAF_only_TR))
kGAF_only_TR<- kGAF_only_TR %>% select("key","SYMBOL",columns_with_AF)
kGAF_only_TR <- kGAF_only_TR%>% select(-c(grep("gnomADe", colnames(kGAF_only_TR))))

# this means it is only in TR but not in SW
# now , rs9332377 only left rs9332377_A, and this has no clinical annotation from PharmaGKB (only rs9332377_T has.)
kGAF_only_TR  <- kGAF_only_TR%>% filter((is.na(ALT_AF_SW)))


# adding drug annotation
kGAF_only_TR_with_drug <- merge(kGAF_only_TR,Target_variants,by = "key")
# rs9332377_A this is shown as both rs9332377_A and rs9332377_T in TR

#write.table(kGAF_only_TR,file = "kGAF_only_TR.xlsx",sep = "\t",quote = FALSE)
write.xlsx(kGAF_only_TR_with_drug, "kGAF_only_TR.xlsx")

# Now remove the rare variants in TR and SW and only keep shared variants 
kGAF <- kGAF %>% filter((!(is.na(ALT_AF_TR) | is.na(ALT_AF_SW))))

rownames(kGAF)<- kGAF$key



#########create a copy of kGAF and plottting
kGAF_test <- kGAF

# Selecting columns by name that contain "AF"
columns_with_AF <- grep("AF", names(kGAF_test))
kGAF_test<- kGAF_test %>% select(columns_with_AF )
kGAF_test <- kGAF_test %>% select(-c(grep("gnomADe", colnames(kGAF_test))))
### TR contains a very rare variants, not present in 1000G

kGAF_test[kGAF_test == "-"] <- 0
kGAF_test <- as.data.frame(kGAF_test)
# Identify columns that contain "-"
kGAF_test
#rows_to_exclude <- apply(kGAF, 1, function(x) any(x == "-"))

# Exclude these columns from kGAF
#kGAF <- kGAF[!rows_to_exclude,]

kGAF_test[] <- lapply(kGAF_test, function(x) as.numeric(as.character(x)))

#kGAF_test <- t(kGAF_test)
pheatmap(kGAF_test,cluster_cols = TRUE,cluster_rows = TRUE,border_color = NA,annotation_row = annotation_row)

write.table(kGAF_test, file= 'kGAF_test.txt',quote = FALSE,sep = '\t', row.names = FALSE)

