# Package loading ------------
library(data.table)
library(dplyr)
library(tidytext)
library(R.utils)
library(ggsci)
library(ggplot2)
dev.off()
library(patchwork)

# Reading files =============
setwd('/Users/xiyas/pgsc_calc/vcf_turkish-with-merge')
var_info <- fread('function_biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz',header = FALSE,sep=' ')
var_info <-distinct(var_info)
colnames(var_info) <- c('CHROM','POS','REF','ALT','rs_id','Consequence','type','IMPACT','MAX_AF','MAX_AF_POPS')
var_info$Consequence <- sapply(strsplit(var_info$Consequence,"&"), `[`, 1)
var_info <- var_info %>%
  mutate(rs_id_category = case_when(
    rs_id == "." ~ "Novel",
    TRUE ~ "Known"
  ))

var_info <- var_info %>%
  mutate(MAX_AF_Category = case_when(
    MAX_AF == "." ~ "No public MAX_AF",
    MAX_AF < 0.01 ~ "Public MAX_AF < 0.01",
    (0.01 < MAX_AF) & (MAX_AF < 0.05) ~ "0.01 < Public MAX_AF < 0.05",
    MAX_AF >= 0.05 ~ "Public MAX_AF >= 0.05"
  ))


#custom_colors <- c("No_public AF" = "orange2", "Public AF < 0.05" = "chartreuse4", "Public AF >= 0.05" = "royalblue3")
custom_colors <- c("deletion" = "tomato1", "SNV" = "chartreuse4", "insertion" = "royalblue3")
legend_labels <- c("insertion" = "INS", "deletion" = "DEL", 'SNV' = 'SNV')

## Get a copy dataframe ============
var_info_test <- var_info
# Reorder the levels of 'Consequence' based on the count in descending order
# Create a summary of the data, ordering it by the count of 'Consequence' in descending order
consequence_counts <- var_info_test %>%
  group_by(rs_id_category, Consequence,type) %>%
  summarize(Count = n())
var_info_test$Consequence <- factor(var_info_test$Consequence, levels =consequence_counts$Consequence)

########Frequency plot. Figure 4A ===============
freq <- fread('freq_biallelic.rm.missing.frq.gz')
colnames(freq) <- c('CHROM','POS','N_ALLELS','N_CHR','REF_AF','ALT_AF')
freq <- freq %>%
  mutate(
    REF = sapply(strsplit(REF_AF, ":", fixed = TRUE), `[`, 1),
    REF_AF = as.numeric(sapply(strsplit(REF_AF, ":", fixed = TRUE), `[`, 2)),
    ALT = sapply(strsplit(ALT_AF, ":", fixed = TRUE), `[`, 1),
    ALT_AF = as.numeric(sapply(strsplit(ALT_AF, ":", fixed = TRUE), `[`, 2))
  )

merge_table_ori <-merge(freq,var_info,all=TRUE,by = c('CHROM','POS','REF','ALT'))


# Define the conditions and group labels
merge_table_ori <- merge_table_ori %>%
  mutate(group = case_when(
    ALT_AF %in% c(0.00181818) ~ "AC =1",
    ALT_AF %in% c(0.00363636) ~ "AC=2",
    ALT_AF %in% c(0.00545455) ~ "AC=3",
    ALT_AF %in% c(0.00727273) ~ "AC=4",
    ALT_AF %in% c(0.00909091) ~ "AC=5",
    (ALT_AF > 0.01 & ALT_AF <= 0.05) ~ "0.01<AF<=0.05",
    (ALT_AF > 0.05 & ALT_AF <= 0.15) ~ "0.05<AF<=0.15",
    (ALT_AF > 0.15 & ALT_AF <= 0.25) ~ "0.15<=AF<=0.25",
    (ALT_AF > 0.25 & ALT_AF <= 0.35) ~ "0.25<=AF<=0.35",
    (ALT_AF > 0.35 & ALT_AF <= 0.45) ~ "0.35<=AF<=0.45",
    (ALT_AF > 0.45 & ALT_AF <= 0.55) ~ "0.45<=AF<=0.55",
    (ALT_AF > 0.55 & ALT_AF <= 0.65) ~ "0.55<=AF<=0.65",
    (ALT_AF > 0.65 & ALT_AF <= 0.75) ~ "0.65<=AF<=0.75",
    (ALT_AF > 0.75 & ALT_AF <= 0.85) ~ "0.75<=AF<=0.85",
    (ALT_AF > 0.85 & ALT_AF <= 0.95) ~ "0.85<=AF<=0.95",
    (ALT_AF > 0.95 & ALT_AF < 1) ~ "0.95<AF<1",
    ALT_AF %in% c(1) ~ "AF=1",
    TRUE ~ "Other"
  )) 


desired_levels <- c(
  "AC =1", "AC=2", "AC=3","AC=4","AC=5",
  "0.01<AF<=0.05", "0.05<AF<=0.15", "0.15<=AF<=0.25",
  "0.25<=AF<=0.35", "0.35<=AF<=0.45", "0.45<=AF<=0.55", "0.55<=AF<=0.65",
  "0.65<=AF<=0.75", "0.75<=AF<=0.85", "0.85<=AF<=0.95", "0.95<AF<1", "AF=1","Other"
)
# Reorder the "group" variable using the desired levels
merge_table_ori$group <- factor(merge_table_ori$group, levels = desired_levels)

# get merge_table =================
### remove chrX and chrY , because the allele frequency calculation weird
merge_table <- merge_table_ori[merge_table_ori$CHROM != 'chrX' & merge_table_ori$CHROM != 'chrY' ,]
merge_table$MAX_AF_Category <- factor(merge_table$MAX_AF_Category, levels = c(
  "Public MAX_AF >= 0.05","0.01 < Public MAX_AF < 0.05", "Public MAX_AF < 0.01", "No public MAX_AF"
))
# Modify specific legend labels
#legend_labels <- c("insertion" = "INS", "deletion" = "DEL", 'SNV' = 'SNV')
write.table(merge_table, file= 'merge_table.txt',quote = FALSE,sep = '\t', row.names = FALSE)

# Figure 4A ggplot ============
ggplot(merge_table, aes(x = group,fill=MAX_AF_Category)) +
  geom_bar(aes(y = ..count.. / 1000000)) +
  labs(x = "Allele Frequency of Turkish cohort", y = "No.of variants (10^6)") +
  scale_fill_jco()+
  #scale_fill_manual(labels = legend_labels) +
  facet_wrap(~ rs_id_category) +
  theme_minimal()+
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # Adjust font sizes
    axis.text.x = element_text(size = 16),     # Font size for x-axis labels
    axis.title.x = element_text(size = 18),    # Font size for x-axis title
    axis.text.y = element_text(size = 16),     # Font size for y-axis labels
    axis.title.y = element_text(size = 18),    # Font size for y-axis title
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Allele_frequencies_Figure4A.pdf",height = 10,width =20)

# Figure 4B ===========
var_info_test <- var_info_test %>%
  mutate(conse_group = case_when(
    Consequence == 'intron_variant' ~ "intronic",
    Consequence == 'intergenic_variant' ~ "intergenic",
    Consequence %in% c('upstream_gene_variant', 'downstream_gene_variant') ~ "regulatory",
    Consequence == '3_prime_utr_variant' ~ "3'-UTR",
    Consequence == '5_prime_utr_variant' ~ "5'-UTR",
    Consequence == 'frameshift_variant ' ~ "frameshift-variant ",
    Consequence == 'missense_variant' ~ "missense",
    Consequence == 'synonymous_variant' ~ "synonymous",
    grepl('splice', Consequence) ~ "splice-region",  # Checking if 'splice' is in Consequence
    Consequence == 'stop_gained' ~ "stop-gain",
    Consequence == 'start_lost' ~ "start-loss",
    Consequence == 'stop_lost' ~ "stop-loss",
    TRUE ~ "Other"
  ))
# Calculate the counts of each Consequence by 'type'
consequence_counts <- var_info_test %>%
  group_by(conse_group) %>%
  summarize(Count = n())%>%
  arrange(desc(Count))

order <- var_info_test %>%
  group_by(conse_group) %>%
  summarize(Count = n())%>%
  arrange(desc(Count))
# Calculate the proportion of 'Novel' within each Consequence
consequence_proportions <- var_info_test %>%
  group_by(conse_group, type, rs_id_category) %>%
  summarize(Count = n()) %>%
  filter(rs_id_category == 'Novel') %>%
  arrange(desc(Count))

MAX_AF_0.05_proportions <- var_info_test %>%
  filter(MAX_AF_Category != 'Public MAX_AF >= 0.05') %>%
  group_by(conse_group) %>%
  summarize(Count = n()) %>%
  arrange(desc(Count))

re_consequence_proportions <- left_join(consequence_proportions, order, by = "conse_group")
re_consequence_proportions <- re_consequence_proportion %>%
  mutate(Novel_Ratio = Count.x / Count.y)

re_MAX_proportions <- left_join(MAX_AF_0.05_proportions, order, by = "conse_group")
re_MAX_proportions <- re_MAX_proportions %>%
  mutate(MAX_Ratio = Count.x / Count.y)
# Arrange consequence_proportions to follow the same order as consequence_counts
consequence_counts$conse_group <- factor(consequence_counts$conse_group, levels =unique(order$conse_group))
re_consequence_proportions$conse_group <- factor(re_consequence_proportions$conse_group, levels =unique(order$conse_group))
re_MAX_proportions$conse_group <- factor(re_MAX_proportions$conse_group, levels =unique(order$conse_group))

# Figure 4B ggplot ============
# Create the facet grid plot with two separate facets for counts and Novel ratios
gg_1 <- ggplot(consequence_counts, aes(x = conse_group, y = Count)) +
  geom_bar(stat = 'identity',fill= "mediumpurple3") +
  #geom_text(aes(label = order$Count), vjust = -0.5, size = 3) + 
  #facet_grid(. ~ Consequence) +
  scale_y_log10() +  # Use a log scale for the y-axis
  #scale_y_continuous(limits = c(0, 5))+
  scale_fill_manual(values = custom_colors,labels = legend_labels) +
  labs(y = "No.of variants (log10)")+
  theme_bw()+ 
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # Adjust font sizes
    axis.text.x = element_blank(),     # Font size for x-axis labels
    axis.title.x = element_blank(),    # Font size for x-axis title
    axis.text.y = element_text(size = 14),     # Font size for y-axis labels
    axis.title.y = element_text(size = 14),    # Font size for y-axis title
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14))

gg_1
# Add a second facet grid for Novel ratios

gg_2 <- ggplot(re_consequence_proportions, aes(x = conse_group, y = Novel_Ratio)) +
  geom_bar(stat = 'identity',fill= "steelblue1") +
  theme_bw()+ 
  #scale_fill_manual(values = "powderblue",labels = legend_labels) +
  labs(x = "Consequences",y= 'Novel (Ratio)')+
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # Adjust font sizes
    axis.text.x = element_text(size = 14),     # Font size for x-axis labels
    axis.title.x = element_blank(),    # Font size for x-axis title
    axis.text.y = element_text(size = 14),     # Font size for y-axis labels
    axis.title.y = element_text(size = 14),    # Font size for y-axis title
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y = element_text(angle = 45, hjust = 1)) 

gg_2

gg_3 <- ggplot(re_MAX_proportions, aes(x = conse_group, y = MAX_Ratio)) +
  geom_bar(stat = 'identity',fill= "orange2") +
  theme_bw()+ 
  #scale_fill_manual(values = "powderblue",labels = legend_labels) +
  labs(x = "Consequences",y='Low Freq (Ratio)')+
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # Adjust font sizes
    axis.text.x = element_blank(),     # Font size for x-axis labels
    axis.title.x = element_blank(),    # Font size for x-axis title
    axis.text.y = element_text(size = 14),     # Font size for y-axis labels
    axis.title.y = element_text(size = 14),    # Font size for y-axis title
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14))+
  theme(axis.title.y = element_text(angle = 45, hjust = 1)) 

gg_3

combined_plot <- gg_1 + gg_3+ gg_2 + plot_layout(nrow = 3, heights = c(10,2,1))  # Adjust the height ratio as needed

# Print the combined plot
print(combined_plot)

# Save the plot as a PDF file
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Figure_4B_Novel_Low_Total.pdf", width = 9, height = 7)

