# Function to read files and process data -----------------
read_and_process_files <- function(path, population) {
  setwd(path)
  files <- list.files(path = path, pattern = "_4_nodup.txt")
  list_temp <- lapply(files, function(file) {
    temp <- read.delim(file)
    # [per-sample] level unique varaints
    temp <- temp[!duplicated(temp$SZAID), ]
    temp$patientID <- rep(file)
    colnames(temp)[26] <- "Sample"
    temp$Population <- population
    return(temp)
  })
  return(do.call(rbind, list_temp))
}

# Function to keep the first pathogenicity data,add condition judgement, rm chrM etc -----------------
customize_temp_data <- function(data) {
  cols <- c("X.CHROM", "POS", "Genotype")
  data$POS <- gsub(" ", "", data$POS)
  data$Variant_info <- apply(data[, cols], 1, paste, collapse = "-")
  data$ClinVar_CLNSIG <- sapply(strsplit(data$ClinVar_CLNSIG, "&"), `[`, 1)
  data$ClinVar_CLNSIG <- sapply(strsplit(data$ClinVar_CLNSIG, "/"), `[`, 1)
  data$Consequence <- sapply(strsplit(data$Consequence, "&"), `[`, 1)
  data <- data %>% filter(X.CHROM != 'chrM')
  data <- data %>%
    mutate(condition = case_when(Zygosity == "Homozygous" ~ "Positive",
                                 Genes %in% AD_genes & Zygosity == "Heterozygous" ~ "Positive",
                                 Genes %in% AR_genes & Zygosity == "Heterozygous" ~ "Carrier",
                                 TRUE ~ "Unsure"))
  return(data)
}

# Function to get the pathogenic variants that suits for writing, to differ with reporting files ------------------
get_P_LP_LoFs <- function(data){
  data <- data %>% filter(X.CHROM != 'chrM') %>% filter(MAX_AF < 0.05 | is.na(MAX_AF)) %>%
  filter(ClinVar_CLNSIG == 'Pathogenic' | ClinVar_CLNSIG == 'Likely_pathogenic' | IMPACT == "HIGH")
  return(data)
}

# Function to get all ClinVar variants -----------------
get_clinvar_variants <- function(data) {
  data <- data %>% filter(grepl("ClinP_LP_var", Database_type))
  return(data)
}
# Function to get other in cillico pLoFs -------------
get_other_sillico_pLoFs <- function(data){
  data <- data %>% filter(grepl('Novel', SZAID))
  return(data)
}

# Function to process ACMG 78 genes findings-----------------
get_ACMG_findings <- function(data) {
  data <- get_P_LP_LoFs(data)
  data <- data %>% filter(final_target_group == 'Basic (for healthy subjects),Usually used for:Health predipositions/Disease risk')
  return(data)
}

# [Cohort-level] unique variants ------------------
get_unique_SZAID <- function(data) {
  data <- data[!duplicated(data$SZAID), ]
  return(data)
}

get_unique_variants <- function(data) {
  data <- data[!duplicated(data$Variant_info), ]
  return(data)
}

# Process Turkish data -----------------
df_temp_turkish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_output_turkish_275/nodup_4_file", "Turkish")
df_temp_turkish <- customize_temp_data(df_temp_turkish)

clinvar_TR <- get_clinvar_variants(df_temp_turkish)
clinvar_TR_unique <- get_unique_SZAID(clinvar_TR)

sillico_pLoFs_TR <- get_other_sillico_pLoFs(df_temp_turkish)
sillico_pLoFs_unique_TR <- get_unique_variants(sillico_pLoFs)

ACMG_TR <- get_ACMG_findings(df_temp_turkish)
ACMG_TR_unique <- get_unique_SZAID(ACMG_TR)



# Process Swedish data -----------------
df_temp_swedish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_output_swedish_101/nodup4_file", "Swedish")
df_temp_swedish <- customize_temp_data(df_temp_swedish)

clinvar_SW <- get_clinvar_variants(df_temp_swedish)
clinvar_SW_unique <- get_unique_SZAID(clinvar_SW)

sillico_pLoFs_SW <- get_other_sillico_pLoFs(df_temp_swedish)
sillico_pLoFs_unique_SW <- get_unique_variants(sillico_pLoFs_SW)

ACMG_SW <- get_ACMG_findings(df_temp_swedish)
ACMG_SW_unique <- get_unique_SZAID(ACMG_SW)



# Combine the two data frames into one -----------------
#combined_df <- rbind(df_clinvar_turkish, df_clinvar_swedish)
#write.table(combined_df,file = "combined_df_clinvar.txt",quote = FALSE,col.names = FALSE,sep = "\t")
#combined_df$Population <- factor(combined_df$Population, levels = c("Turkish", "Swedish"))
#clinvar_genes_count_combined <- rbind(clinvar_genes_count_turkish, clinvar_genes_count_swedish)
