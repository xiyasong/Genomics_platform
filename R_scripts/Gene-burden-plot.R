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

df_temp_turkish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_results_TR_275_2025", "Turkish") %>% 
  customize_temp_data()
# ClinVar variants
clinvar_TR <- get_clinvar_variants(df_temp_turkish)
clinvar_TR_unique <- clinvar_TR %>% get_unique_variants()  # 667 × 163

df_temp_swedish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_results_SW_101_2025", "Swedish") %>% 
  customize_temp_data()

# (rare variant bining set methods 1: TR)
# Function to process data based on selected method and cohort
process_variants <- function(data, cohort, method) {
  # Select cohort data
  df <- if(cohort == "TR") df_temp_turkish else df_temp_swedish
  
  # Apply filtering based on selected method
  temp_table <- switch(method,
    # Way 1: ClinVar variants (non-conflicting)
    {
      get_clinvar_variants(df) %>% 
        filter(
          MAX_AF <= 0.05,
          ClinVar_CLNSIG != "Conflicting_classifications_of_pathogenicity"
        )
    },
    # Way 2: High impact variants
    {
      df %>% 
        filter(
          MAX_AF <= 0.05,
          IMPACT == 'HIGH'
        )
    },
    # Way 3: Variants with reason specified
    {
      df %>% 
        filter(
          MAX_AF <= 0.05,
          Reason != ""
          #IMPACT == 'HIGH'|am_class == 'likely_pathogenic'
        )
    }
  ) %>%  # This pipe applies to ALL methods
    add_count(Variant_info, name = "Freq") %>%
    distinct(Variant_info, .keep_all = TRUE) %>%
    select(Freq, everything())
  
  return(temp_table)
}
# Usage examples by setting three strategies ----------------------
# For Turkish cohort, way 1
temp_table <- process_variants(cohort = 'SW', method = 3)
temp_table <- temp_table %>%
  mutate(group = case_when(
    Freq == 1 ~ "count(1)",
    Freq > 1 & Freq <= 4 ~ "count(2-4)",
    Freq > 4 & Freq <= 7 ~ "count(5-7)", 
    Freq > 7 & Freq <= 10 ~ "count(8-10)",
    Freq > 10 & Freq <= 15 ~ "count(11-15)",
    Freq > 15 & Freq <= 20 ~ "count(16-20)",
    Freq > 20 & Freq <= 30 ~ "count(21-30)",
    Freq > 30 & Freq <= 40  ~ "count(>30)",
    Freq > 40 ~ "count(>40)"
  )) %>%  
  arrange(-desc(Genes), -desc(group)) %>% 
  as_tibble()

gene_counts_TR<- temp_table%>%
  group_by(Genes) %>%
  summarise(Total_Freq = sum(Freq)) %>%
  arrange(desc(Total_Freq)) 

gene_selected_TR = gene_counts_TR$Genes[1:30]

desired_levels <- c(
 "count(1)", "count(2-4)", "count(5-7)", "count(8-10)",
 "count(11-15)", "count(16-20)", "count(21-30)", "count(>30)","count(>40)"
)

# Get all 10 colors from RdYlBu palette
all_colors <-  RColorBrewer::brewer.pal(n = 10, name = "RdYlBu") 
# Reorder colors to match your specified mapping:
# count(1)=6, count(2~5)=7, count(5~10)=8, count(10~15)=9, count(15~20)=10,
# count(20~25)=5, count(25~30)=4, count(30~35)=3, count(35~40)=2, count(>40)=1
ordered_colors <- all_colors[c(6,7,9,10,5,4,3,2,1)]
plot_data_TR <- temp_table %>% 
  dplyr::select(Genes, group) %>% 
  filter(Genes %in% gene_selected_TR) %>% 
  group_by(Genes, group) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  mutate(group = factor(group, levels = desired_levels))

# Create confidence color settings
confidence_setting <- ifelse(gene_selected_TR %in% High_genes$Genes, "darkblue", "black")

max_total_count <- plot_data_TR %>% 
  group_by(Genes) %>% 
  summarise(total_count = sum(count)) %>% 
  pull(total_count) %>% 
  max() %>% 
  ceiling() +5  # Rounds up to nearest integer and adds 2


# Generate the plot
plot_gene_burden <- plot_data_TR %>% 
  ggplot(aes(x = Genes, y = count, fill = group)) +
  geom_bar(stat = "identity", color = "white", alpha = 0.9) +
  scale_x_discrete(limits = gene_selected_TR) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, max_total_count)) +
  # selection 1
  scale_fill_manual(values = ordered_colors,
    breaks = desired_levels) +
  #scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,8,9,10,5,4,3,1)]) +
  ggtitle("SW Carrier frequency of pLoF+p-risk variants") +
  ylab("No. unique variants") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5),
            size = 3.5) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_line(size = 0.5, colour = "grey"),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, 
                             colour = confidence_setting, 
                             angle = 315, 
                             vjust = 0.5, 
                             hjust = 0)
  )
plot_gene_burden

```
## 7.1 save plot
```{r save}

ggsave(plot_gene_burden, 
       filename = "/Users/xiyas/V2_Genome_reporting/Plots/gene-burden-3-SW.pdf",
       width = 680/72,
       height = 360/72)
