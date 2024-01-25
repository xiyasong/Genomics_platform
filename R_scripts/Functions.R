# The functions and main to generate these dataframes:
# df_clinvar = clinvar_TR
# test = clinvar_TR_unique
# df_puta =sillico_pLoFs_TR
# combined_df = combined_df
# sorted_tab_old = clinvar_TR_unique + freq calculation


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
  #cols <- c("X.CHROM", "POS", "Genotype")
  cols <- c("X.CHROM", "POS", "REF",'ALT')
  #sillico_pLoFs_TR$SZAID_assign <- apply(sillico_pLoFs_TR[, cols], 1, paste, collapse = "-")
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
  data <- data %>%
    mutate(MAX_AF_Category = case_when(
      is.na(MAX_AF) ~ "No public MAX_AF",
      MAX_AF < 0.01 ~ "Public MAX_AF < 0.01",
      (0.01 < MAX_AF) & (MAX_AF < 0.05) ~ "0.01 < Public MAX_AF < 0.05",
      MAX_AF >= 0.05 ~ "Public MAX_AF >= 0.05"
    ))
  return(data)
}

# Function to get the pathogenic variants that suits for writing, to differ with reporting files ------------------
get_P_LP_LoFs <- function(data){
  data <- data %>% filter(X.CHROM != 'chrM') %>% filter(MAX_AF < 0.05 | is.na(MAX_AF)) %>%
    filter(ClinVar_CLNSIG == 'Pathogenic' | ClinVar_CLNSIG == 'Likely_pathogenic' | IMPACT == "HIGH")
  return(data)
}

get_P_LP <- function(data){
  data <- data %>% filter(X.CHROM != 'chrM') %>% filter(MAX_AF < 0.05 | is.na(MAX_AF)) %>%
    filter(ClinVar_CLNSIG == 'Pathogenic' | ClinVar_CLNSIG == 'Likely_pathogenic')
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

# for novel pLoFs that do not have SZAID --------------
get_unique_variants <- function(data) {
  data <- data[!duplicated(data$Variant_info), ]
  return(data)
}

#  functions to get sorted_tab --------------
library(dplyr)

get_sorted_tab <- function(data, unique_data) {
  sorted_tab <- data %>%
    group_by(SZAID, Zygosity) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    merge(unique_data, by = 'SZAID', all.x = TRUE)
  
  sorted_tab_old <- data %>%
    group_by(SZAID) %>%
    summarise(Freq = n()) %>%
    arrange(desc(Freq)) %>%
    merge(unique_data, by = 'SZAID', all.x = TRUE) 
  
  return(list(sorted_tab = sorted_tab, sorted_tab_old = sorted_tab_old))
}

get_sorted_tab_sillico <- function(data, unique_data) {
  sorted_tab_sillico <- data %>%
    group_by(Variant_info, Zygosity) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    merge(unique_data, by = 'Variant_info', all.x = TRUE)
  
  sorted_tab_old_sillico <- data %>%
    group_by(Variant_info) %>%
    summarise(Freq = n()) %>%
    arrange(desc(Freq)) %>%
    merge(unique_data, by = 'Variant_info', all.x = TRUE) 
  
  return(list(sorted_tab_sillico = sorted_tab_sillico, sorted_tab_old_sillico = sorted_tab_old_sillico))
}