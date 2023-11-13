# Function to read files and process data -----------------
read_and_process_files <- function(path, population) {
  setwd(path)
  files <- list.files(path = path, pattern = "_4_nodup.txt")
  
  list_temp <- lapply(files, function(file) {
    temp <- read.delim(file)
    cols <- c("X.CHROM", "POS", "Genotype")
    temp$POS <- gsub(" ", "", temp$POS)
    temp$Variant_info <- apply(temp[, cols], 1, paste, collapse = "-")
    temp <- temp[!duplicated(temp$SZAID), ]
    temp$Variant_info <- apply(temp[, cols], 1, paste, collapse = "-")
    temp$patientID <- rep(file)
    colnames(temp)[26] <- "Sample"
    return(temp)
  })
  
  return(do.call(rbind, list_temp))
}

# Function to normalize pathogenicity data and add condition judgement -----------------
customize_pathogenicity <- function(data) {
  data$ClinVar_CLNSIG <- sapply(strsplit(data$ClinVar_CLNSIG, "&"), `[`, 1)
  data$ClinVar_CLNSIG <- sapply(strsplit(data$ClinVar_CLNSIG, "/"), `[`, 1)
  data <- data %>% filter(X.CHROM != 'chrM')
  data <- data %>%
    mutate(condition = case_when(Zygosity == "Homozygous" ~ "Positive",
                                 Genes %in% AD_genes & Zygosity == "Heterozygous" ~ "Positive",
                                 Genes %in% AR_genes & Zygosity == "Heterozygous" ~ "Carrier",
                                 TRUE ~ "Unsure"))
  return(data)
}

# Function to process ACMG 78 genes findings-----------------
process_ACMG_findings <- function(data) {
  df_temp_ACMG <- data %>% filter(final_target_group == 'Basic (for healthy subjects),Usually used for:Health predipositions/Disease risk')
  df_temp_ACMG_clin_lof <- normalize_pathogenicity(df_temp_ACMG) %>%
    filter(ClinVar_CLNSIG == 'Pathogenic' | ClinVar_CLNSIG == 'Likely_pathogenic' | IMPACT == "HIGH") %>%
    filter(X.CHROM != 'chrM') %>%
    filter(MAX_AF < 0.05 | is.na(MAX_AF))
  df_temp_ACMG_clin_lof_unique <- df_temp_ACMG_clin_lof[!duplicated(df_temp_ACMG_clin_lof$SZAID), ]
  ACMG_plofs <- df_temp_ACMG_clin_lof_unique %>% filter(grepl('Novel', SZAID))
  ACMG_plofs$Consequence <- sapply(strsplit(ACMG_plofs$Consequence, "&"), `[`, 1)
  # Return as a list
  return(list(df_temp_ACMG = df_temp_ACMG,
              df_temp_ACMG_clin_lof = df_temp_ACMG_clin_lof, 
              df_temp_ACMG_clin_lof_unique = df_temp_ACMG_clin_lof_unique,
              ACMG_plofs = ACMG_plofs))
}

# Function to process ClinVar variants -----------------
process_clinvar_variants <- function(data, population) {
  df_clinvar <- data %>% filter(grepl("ClinP_LP_var", Database_type))
  df_clinvar <- df_clinvar[!duplicated(df_clinvar$SZAID), ]
  clinvar_genes_count <- df_clinvar %>%
    group_by(Genes, ClinVar_CLNSIG) %>%
    summarize(Count = n()) %>%
    arrange(desc(Count))
  df_clinvar$Population <- population
  return(list(df_clinvar, clinvar_genes_count))
}

# Process Turkish data -----------------
df_temp_turkish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_output_turkish_275/nodup_4_file", "Turkish")
df_temp_turkish <- customize_pathogenicity(df_temp_turkish)
ACMG_results_turkish <- process_ACMG_findings(df_temp_turkish)
clinvar_turkish <- process_clinvar_variants(df_temp_turkish, "Turkish")
# Access the components using [[]]
df_clinvar_turkish <- clinvar_turkish[[1]]
clinvar_genes_count_turkish <- clinvar_turkish[[2]]

# Process Swedish data -----------------
df_temp_swedish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_output_swedish_101/nodup4_file", "Swedish")
df_temp_swedish <- customize_pathogenicity(df_temp_swedish)
clinvar_swedish <- process_clinvar_variants(df_temp_swedish, "Swedish")
# Access the components using [[]]
df_clinvar_swedish <- clinvar_swedish[[1]]
clinvar_genes_count_swedish <- clinvar_swedish[[2]]

# Combine the two data frames into one -----------------
combined_df <- rbind(df_clinvar_turkish, df_clinvar_swedish)
combined_df$Population <- factor(combined_df$Population, levels = c("Turkish", "Swedish"))
#clinvar_genes_count_combined <- rbind(clinvar_genes_count_turkish, clinvar_genes_count_swedish)
