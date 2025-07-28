# The functions and main to generate these dataframes:
# df_clinvar = clinvar_TR
# test = clinvar_TR_unique
# df_puta =sillico_pLoFs_TR
# combined_df = combined_df
# sorted_tab_old = clinvar_TR_unique + freq calculation


# Function to read files and process data -----------------
read_and_process_files <- function(path, population) {
  # Set working directory (optional, but ensure path is correct)
  if (!dir.exists(path)) stop("Directory does not exist: ", path)
  # List all files ending with '_4_nodup.txt'
  files <- list.files(path = path, pattern = "_4_nodup\\.txt$", full.names = TRUE)
  if (length(files) == 0) stop("No files matching '_4_nodup.txt' found in: ", path)
  # Read each file and standardize columns
  list_temp <- lapply(files, function(file) {
    temp <- read.delim(file, stringsAsFactors = FALSE)
    # Ensure 'SZAID' exists (avoid errors if column is missing)
    if (!"SZAID" %in% colnames(temp)) {
      warning("File '", basename(file), "' is missing column 'SZAID'. Skipping deduplication.")
    } else {
      temp <- temp[!duplicated(temp$SZAID), ]  # Deduplicate by SZAID
    }
    
    # Add metadata columns
    temp$patientID <- basename(file)  # Use basename() to avoid full paths
    temp$Population <- population
    
    # Rename column 26 to 'Sample' (if it exists)
    if (ncol(temp) >= 26) {
      colnames(temp)[26] <- "Sample"
    } else {
      warning("File '", basename(file), "' has fewer than 26 columns. 'Sample' not assigned.")
    }
    
    return(temp)
  })
  
  # Standardize columns across all data frames before binding
  all_cols <- unique(unlist(lapply(list_temp, colnames)))
  list_temp_std <- lapply(list_temp, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA  # Fill missing columns with NA
    }
    return(df)
  })
  
  # Combine all data frames
  final_df <- do.call(rbind, list_temp_std)
  return(final_df)
}


# Function to keep the first pathogenicity data,add condition judgement, rm chrM etc -----------------
customize_temp_data <- function(data) {
  # Basic columns to create Variant_info
  cols <- c("X.CHROM", "POS", "REF", "ALT")
  
  # Clean up POS and generate variant ID
  data$POS <- gsub(" ", "", data$POS)
  data$Variant_info <- apply(data[, cols], 1, paste, collapse = "-")
  
  # Clean up annotation columns
  data$ClinVar_CLNSIG <- sapply(strsplit(data$ClinVar_CLNSIG, "&"), `[`, 1)
  data$ClinVar_CLNSIG <- sapply(strsplit(data$ClinVar_CLNSIG, "/"), `[`, 1)
  data$Consequence <- sapply(strsplit(data$Consequence, "&"), `[`, 1)
  
  # Filter out mitochondrial chromosome
  data <- data %>% filter(X.CHROM != 'chrM')
  
  # Determine condition based on Inheritances and Zygosity
  data <- data %>%
    mutate(
      condition = case_when(
        grepl("recessive", Inheritances, ignore.case = TRUE) & Zygosity == "Homozygous" ~ "Positive",
        grepl("recessive", Inheritances, ignore.case = TRUE) & Zygosity == "Heterozygous" ~ "Carrier",
        grepl("dominant", Inheritances, ignore.case = TRUE) & Zygosity == "Heterozygous" ~ "Positive",
        TRUE ~ "Unsure"
      )
    )
  
  # Determine MAX_AF category
  data <- data %>%
    mutate(
      MAX_AF_Category = case_when(
        is.na(MAX_AF) ~ "No public MAX_AF",
        MAX_AF < 0.01 ~ "Public MAX_AF < 0.01",
        (0.01 <= MAX_AF) & (MAX_AF < 0.05) ~ "0.01 <= Public MAX_AF < 0.05",
        MAX_AF >= 0.05 ~ "Public MAX_AF >= 0.05"
      )
    )
  
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

# for novel pLoFs that do not have SZAID --------------
get_unique_variants <- function(data) {
  data <- data[!duplicated(data$Variant_info), ]
  return(data)
}

#  functions to get sorted_tab --------------
library(dplyr)


read_GWAS <- function(path) {
  setwd(path)
  files <- list.files(path = path, pattern = "_GWAS.txt")
  list_temp <- lapply(files, function(file) {
    temp <- read.delim(file)
    temp <- temp[!duplicated(temp$SZAvarID), ]
    temp$patientID <- rep(file)
    return(temp)
  })
  return(do.call(rbind, list_temp))
}
read_pharmaco <- function(path) {
  setwd(path)
  files <- list.files(path = path, pattern = "_pharmaco.txt")
  list_temp <- lapply(files, function(file) {
    temp <- read.delim(file)
    temp <- temp[!duplicated(temp$SZAvarID), ]
    temp$patientID <- rep(file)
    return(temp)
  })
  return(do.call(rbind, list_temp))
}

