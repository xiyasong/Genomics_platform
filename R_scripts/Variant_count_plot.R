setwd("/Users/xiyas/V2_Genome_reporting/python_output_file_swedish/")
library(dplyr)
library(ggsci)
library(ggplot2)
# Set working directory
setwd("/Users/xiyas/V2_Genome_reporting/python_results_TR_275_2025/")

#setwd("/Users/xiyas/V2_Genome_reporting/python_results_SW_101_2025/")
# Get file list
files <- list.files(path = ".", pattern = "_4_nodup.txt")
# Define categories
categories <- c("Health Predisposition", "Newborn screening", "Carrier screening", 
                "Heriditary Cancer screening", "Mono-rare pairs", "Clinvar all pairs")
# Create a function to filter data
filter_variants <- function(data, db_pattern, category_pattern) {
  data %>% 
    filter(grepl(db_pattern, Database_type) | grepl(db_pattern, SZAID)) %>% 
    filter(grepl(category_pattern, final_target_group)) %>% 
    pull(SZAID) %>% 
    unique() %>% 
    length()
}
# Create mapping for category patterns
category_patterns <- c(
  "Health predipositions",
  "Newborn-screening",
  "Carrier-screening",
  "Heriditary-cancer",
  "Expanded_mono_rare_diseases",
  "Expanded_clinvar_mono_diseases"
)

# Initialize results dataframes
results_clinvar <- matrix(0, nrow = length(files), ncol = length(categories))
results_puta <- matrix(0, nrow = length(files), ncol = length(categories))

# Process each file
for (i in seq_along(files)) {
  data <- read.delim(files[i])
  
  # Process each category
  for (j in seq_along(categories)) {
    results_clinvar[i, j] <- filter_variants(data, "ClinP_LP_var", category_patterns[j])
    results_puta[i, j] <- filter_variants(data, "NovelTrans", category_patterns[j])
  }
}
# Convert results to dataframes
df_count_clinvar <- as.data.frame(results_clinvar)
df_count_puta <- as.data.frame(results_puta)

# Set column names
colnames(df_count_clinvar) <- categories
colnames(df_count_puta) <- categories

# Prepare data for plotting
df_count_clinvar <- melt(df_count_clinvar)
df_count_clinvar$group <- "Clinvar P/LP variants"

df_count_puta <- melt(df_count_puta)
df_count_puta$group <- "Putative risk variants"

df_sum <- rbind(df_count_clinvar, df_count_puta)

# Create plots
# Clinvar plot
# For the Clinvar plot
p1 <- ggplot(df_count_clinvar, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_fill_lancet() +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40)) +
  ylab("Gene panels") +
  xlab("All Clinvar P/LP variants split by gene panels")

# For the Putative variants plot
p2 <- ggplot(df_count_puta, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_fill_lancet() +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40))+
  ylab("Gene panels") +
  xlab("All putative risk variants split by gene panels")

# For the Combined plot
p3 <- ggplot(df_sum, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges() +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  ) +
  scale_fill_lancet() +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40))+
  ylab("Gene panels") +
  xlab("No.of variants on each TR sample") +
  facet_wrap(~group, ncol = 1)

p3

ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Variant_count_ridgeline_plot_275_TR_2025.pdf",height = 6,width = 8)

#ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Variant_count_ridgeline_plot_101_SW_2025.pdf",height = 6,width = 8)




  