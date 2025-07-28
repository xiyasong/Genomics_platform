library(dplyr)
library(ggplot2)
library(reshape2)

# Combine Turkish and Swedish pharmaco data, excluding homo_ref: using level of evidence 1 and 2
combined_pharmaco <- bind_rows(df_pharmaco_TR_remove_4, df_pharmaco_SW_remove_4) %>% 
  filter(Level.of.Evidence != 3) %>% 
  mutate(
    # Create unique combination ID: variant + genotype
    variant_genotype = paste(rsID, pharma_genotype, sep = "-")
  ) %>%
  select(variant_genotype, Drug.s., patientID, Population, Level.of.Evidence)
# 计算 TR 的 median
a_TR <- combined_pharmaco %>%
  filter(Population == "Turkish") %>%
  count(patientID)
median_TR <- median(a_TR$n)

# 计算 SW 的 median
a_SW <- combined_pharmaco %>%
  filter(Population == "Swedish") %>%
  count(patientID)
median_SW <- median(a_SW$n)

# 显示结果
median_TR
median_SW

# Calculate total patients per population
total_patients_pharma <- combined_pharmaco %>%
  group_by(Population) %>%
  summarise(total = n_distinct(patientID)) %>%
  ungroup()

# Calculate variant-genotype frequencies
variant_freq_pharma <- combined_pharmaco %>%
  group_by(variant_genotype, Drug.s., Population, Level.of.Evidence) %>%
  summarise(n_carriers = n_distinct(patientID), .groups = "drop") %>%
  left_join(total_patients_pharma, by = "Population") %>%
  mutate(freq = n_carriers / total)
# Convert to wide format for difference calculation
freq_wide_pharma <- variant_freq_pharma %>%
  dcast(variant_genotype + Drug.s. + Level.of.Evidence ~ Population, 
        value.var = "freq", fill = 0) %>%
  mutate(
    freq_diff = abs(Turkish - Swedish),
    direction = ifelse(Turkish > Swedish, "Higher in TR", "Higher in SW")
  ) %>%
  arrange(desc(freq_diff))

# Filter out Swedish=0 and get top 20
top_variants_pharma <- freq_wide_pharma %>% 
  filter(Swedish != 0) %>%
  head(20) %>% 
  pull(variant_genotype)

plot_data_pharma <- freq_wide_pharma %>% 
  filter(variant_genotype %in% top_variants_pharma)

# Custom label wrapping function
wrap_labels <- function(x, width = 15) {
  sapply(x, function(s) {
    paste(strwrap(s, width = width), collapse = "\n")
  })
}

# Extract first drug from semicolon-separated list
plot_data_pharma$Drug.s. <- sapply(strsplit(plot_data_pharma$Drug.s., ";"), function(x) trimws(x[1]))

# Create faceted heatmap
pgx_compare <- ggplot(plot_data_pharma, aes(x = wrap_labels(Drug.s.), y = variant_genotype)) +
  geom_tile(aes(fill = freq_diff), color = "white", linewidth = 0.3) +
  geom_point(
    aes(shape = direction),
    size = 3,
    color = "black",
    position = position_nudge(y = -0.2)
  ) +
  facet_grid(Level.of.Evidence ~ ., scales = "free_y", space = "free") +
  scale_fill_gradient(
    name = "Frequency Difference",
    low = "#E6E6FA",  # Light purple
    high = "#4B0082",  # Dark purple
    limits = c(0.09, 0.4)
  ) +
  scale_shape_manual(
    name = "Direction",
    values = c("Higher in TR" = 16, "Higher in SW" = 17),
    guide = guide_legend(override.aes = list(size = 5))  # larger shapes in legend

  ) +
  labs(
    title = "Top 20 PGx-SNPs with Largest Frequency Differences",
    x = "Associated Drug (primary only)",
    y = "Variant-Genotype Combination"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text.y = element_text(angle = 0),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
pgx_compare
ggsave(pgx_compare,
       filename = "/Users/xiyas/V2_Genome_reporting/Plots/pgx_compare.pdf",width = 897/72,height = 530/72)
library(writexl)

write_xlsx(
  list(
    "Traits" = freq_wide,
    "PGx" = freq_wide_pharma
  ),
  path = "/Users/xiyas/V2_Genome_reporting/data/trait_pgx_carrier_freq_compare.xlsx")


########## not comparing differencies, but rather the top carrier finidngs =====================
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(viridis)
library(patchwork)
library(tidytext)  # 用于reorder_within和scale_x_reordered函数

# 首先，确保人群标签一致性
pgx_data <- variant_freq_pharma %>%
  mutate(Population = ifelse(Population == "Turkish", "TR", 
                             ifelse(Population == "Swedish", "SW", Population)))

# 为第一个图表，分别获取每个人群的顶变异
# 获取TR的前10个变异
tr_top_pgx <- pgx_data %>%
  filter(Population == "TR") %>%
  arrange(desc(freq)) %>%
  head(10) %>%
  pull(variant_genotype)

# 获取SW的前10个变异
sw_top_pgx <- pgx_data %>%
  filter(Population == "SW") %>%
  arrange(desc(freq)) %>%
  head(10) %>%
  pull(variant_genotype)

# 合并列表，获取唯一变异
combined_top_pgx <- unique(c(tr_top_pgx, sw_top_pgx))

# 筛选只包含这些顶变异的数据
viz_pgx_data <- pgx_data %>%
  filter(variant_genotype %in% combined_top_pgx) %>%
  # 创建更易读的标签，包含药物和变异信息
  mutate(var_display_name = paste0(variant_genotype, "\n", Drug.s.))

# 创建按人群分面排序的条形图
bar_plot_pgx <- viz_pgx_data %>%
  ggplot(aes(x = reorder_within(var_display_name, freq, Population), y = freq, fill = Population)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("SW" = "#3b9ab2", "TR" = "#e1812c")) +
  scale_x_reordered() +
  facet_wrap(~ Population, scales = "free_y") +
  labs(title = "Top Pharmacogenomics Variants by Frequency",
       x = NULL,
       y = "Carrier Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  coord_flip()

bar_plot_pgx
# 2. 按证据级别和药物分类的热力图
heat_plot_pgx <- viz_pgx_data %>%
  ggplot(aes(x = Population, y = fct_reorder(var_display_name, freq), fill = freq)) +
  geom_tile() +
  geom_text(aes(label = Level.of.Evidence), color = "black", size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "plasma", name = "Frequency") +
  labs(title = "PGx Variant Frequency by Evidence Level",
       x = "Population",
       y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")

# 3. 人群间频率比较的散点图
pgx_comparison <- pgx_data %>%
  select(variant_genotype, Drug.s., Population, freq, Level.of.Evidence) %>%
  pivot_wider(names_from = Population, values_from = freq, values_fill = 0) %>%
  filter(!is.na(TR) & !is.na(SW))

scatter_plot_pgx <- ggplot(pgx_comparison, aes(x = TR, y = SW, color = Level.of.Evidence)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_viridis_d(option = "turbo") +
  xlim(0, max(c(pgx_comparison$TR, pgx_comparison$SW)) * 1.05) +
  ylim(0, max(c(pgx_comparison$TR, pgx_comparison$SW)) * 1.05) +
  labs(title = "PGx Variant Frequency: TR vs SW",
       x = "TR Population Frequency",
       y = "SW Population Frequency",
       color = "Evidence Level") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")
scatter_plot_pgx

# 4. 按药物分类的中位频率气泡图
bubble_plot_pgx <- pgx_data %>%
  group_by(Drug.s., Population) %>%
  summarise(max_freq = max(freq),
            median_freq = median(freq),
            variant_count = n_distinct(variant_genotype),
            .groups = "drop") %>%
  ggplot(aes(x = fct_reorder(Drug.s., median_freq), y = median_freq, 
             fill = Population, size = variant_count)) +
  geom_point(shape = 21, alpha = 0.8, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("SW" = "#3b9ab2", "TR" = "#e1812c")) +
  scale_size_continuous(range = c(3, 12)) +
  labs(title = "PGx Variants by Drug",
       x = NULL,
       y = "Median Carrier Frequency",
       size = "Number of\nVariants",
       fill = "Population") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  coord_flip()

# 5. 按证据级别和人群的分布图
ridgeline_plot_pgx <- ggplot(pgx_data, aes(x = freq, fill = Population)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("SW" = "#3b9ab2", "TR" = "#e1812c")) +
  facet_wrap(~Level.of.Evidence, scales = "free_y") +
  labs(title = "Frequency Distribution by Evidence Level",
       x = "Variant Frequency",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top")

# 6. 额外：药物特定的变异频率比较
drug_specific_plot <- pgx_data %>%
  filter(Drug.s. %in% c("nicotine", "warfarin", "clopidogrel")) %>%  # 选择几个常见药物
  ggplot(aes(x = variant_genotype, y = freq, fill = Population)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("SW" = "#3b9ab2", "TR" = "#e1812c")) +
  facet_wrap(~Drug.s., scales = "free_x") +
  labs(title = "Variant Frequencies for Selected Drugs",
       x = "Variant",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top")
drug_specific_plot
# 使用patchwork组合图表
multi_panel_pgx <- (bar_plot_pgx | heat_plot_pgx) / (scatter_plot_pgx | bubble_plot_pgx)
multi_panel_pgx + plot_annotation(
  title = "Pharmacogenomics Variant Frequency Comparison Between TR and SW Populations",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)

# 查看单独的图表
# bar_plot_pgx     # 顶部变异条形图
# heat_plot_pgx    # 按证据级别的热力图
# scatter_plot_pgx # 人群比较散点图
# bubble_plot_pgx  # 按药物的气泡图
# ridgeline_plot_pgx # 按证据级别的分布图
# drug_specific_plot # 特定药物的变异频率
########## finding TR only PGx variants --------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)

# 整理数据，确保人群标签一致性
pgx_data <- variant_freq_pharma %>%
  mutate(Population = ifelse(Population == "Turkish", "TR", 
                             ifelse(Population == "Swedish", "SW", Population)))

# 找出每个人群中存在的变异
tr_variants <- pgx_data %>%
  filter(Population == "TR") %>%
  distinct(variant_genotype) %>%
  pull(variant_genotype)

sw_variants <- pgx_data %>%
  filter(Population == "SW") %>%
  distinct(variant_genotype) %>%
  pull(variant_genotype)

# 找出只在TR中存在，但在SW中不存在的变异
tr_only_variants <- setdiff(tr_variants, sw_variants)

# 获取这些TR特有变异的详细信息
tr_specific_variants <- pgx_data %>%
  filter(Population == "TR" & variant_genotype %in% tr_only_variants) %>%
  arrange(desc(freq)) %>%
  select(variant_genotype, Drug.s., Level.of.Evidence, n_carriers, total, freq)

# 打印TR特有变异的结果
print(paste("Total variants in TR:", length(tr_variants)))
print(paste("Total variants in SW:", length(sw_variants)))
print(paste("Variants only in TR but not in SW:", length(tr_only_variants)))
print("")
print("TR-specific variants (sorted by frequency):")
print(tr_specific_variants)

# 创建一个更美观的输出表格
tr_specific_summary <- tr_specific_variants %>%
  mutate(
    percent_carriers = round(freq * 100, 2),
    variant_info = paste0(variant_genotype, " (", Drug.s., ")")
  ) %>%
  select(Variant = variant_info, 
         Evidence = Level.of.Evidence,
         Carriers = n_carriers,
         Total = total,
         `Frequency (%)` = percent_carriers)

print("Formatted summary of TR-specific variants:")
print(kable(tr_specific_summary, caption = "Pharmacogenomic variants found only in Turkish population"))

# 创建可视化展示这些TR特有变异
tr_spec_plot <- ggplot(tr_specific_variants, aes(x = reorder(variant_genotype, freq), y = freq)) +
  geom_bar(stat = "identity", fill = "#e1812c", alpha = 0.8) +
  geom_text(aes(label = paste0(Drug.s., " (", Level.of.Evidence, ")")), 
            hjust = -0.05, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Pharmacogenomic Variants Found Only in Turkish Population",
       x = "Variant",
       y = "Frequency in TR Population") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold")) +
  coord_flip()

# 显示图表
print(tr_spec_plot)

# 按药物分类的TR特有变异
tr_by_drug <- tr_specific_variants %>%
  group_by(Drug.s.) %>%
  summarise(
    n_variants = n(),
    max_freq = max(freq),
    total_carriers = sum(n_carriers),
    evidence_levels = paste(unique(Level.of.Evidence), collapse = ", ")
  ) %>%
  arrange(desc(max_freq))

print("TR-specific variants grouped by drug:")
print(tr_by_drug)

# 统计文字摘要
cat("\n=== SUMMARY ===\n")
cat(paste("- Turkish population has", length(tr_only_variants), "unique PGx variants not found in Swedish population\n"))
cat(paste("- These TR-specific variants are associated with", length(unique(tr_specific_variants$Drug.s.)), "different drugs\n"))
cat(paste("- The most common TR-specific variant is", tr_specific_variants$variant_genotype[1], 
          "with frequency", round(tr_specific_variants$freq[1] * 100, 2), "% (", tr_specific_variants$Drug.s.[1], ")\n"))
cat(paste("- Evidence levels represented:", paste(unique(tr_specific_variants$Level.of.Evidence), collapse = ", "), "\n"))

