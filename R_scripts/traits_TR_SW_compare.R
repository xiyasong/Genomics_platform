library(dplyr)
library(stringr)

##########Figure 8A --------------------------
combined_traits <- bind_rows(df_traits_TR, df_traits_SW) %>% filter(Zygosity != "homo_ref") %>% 
  mutate(
    # Create unique combination ID
    variant_genotype = paste(variants, Patient.genotypes, sep = "-")
  ) %>%
  select(variant_genotype, Traits.name,patientID, Population,category)

total_patients <- combined_traits %>%
  group_by(Population) %>%
  summarise(total = n_distinct(patientID)) %>%
  ungroup()

variant_freq <- combined_traits %>%
  group_by(variant_genotype,Traits.name, Population,category) %>%
  summarise(n_carriers = n_distinct(patientID), .groups = "drop") %>%
  left_join(total_patients, by = "Population") %>%
  mutate(freq = n_carriers / total)

head(variant_freq)
library(ggplot2)
library(reshape2)

freq_wide <- variant_freq %>%
  dcast(variant_genotype + Traits.name + category ~ Population, value.var = "freq", fill = 0) %>%
  mutate(
    freq_diff = abs(Turkish - Swedish),
    direction = ifelse(Turkish > Swedish, "Higher in TR", "Higher in SW")
  ) %>%
  arrange(desc(freq_diff))

#
top_variants <- head(freq_wide$variant_genotype, 20)
plot_data <- freq_wide %>% filter(variant_genotype %in% top_variants)

# 自定义标签换行函数
wrap_labels <- function(x, width = 15) {
  sapply(x, function(s) {
    paste(strwrap(s, width = width), collapse = "\n")
  })
}

# 绘制热图
trait_compare <- ggplot(plot_data, aes(x = wrap_labels(Traits.name), y = variant_genotype)) +
  geom_tile(aes(fill = freq_diff), color = "white", linewidth = 0.3) +
  geom_point(
    aes(shape = direction),
    size = 3,
    color = "black",
    position = position_nudge(y = -0.2)
  ) +
  facet_grid(category ~ ., 
             scales = "free_y", 
             space = "free",
             labeller = label_wrap_gen(width = 10)) +
  scale_fill_gradient(
    name = "Frequency Difference",
    low = "#E6E6FA",  # 浅紫色
    high = "#4B0082",  # 深紫色
    limits = c(0.2, 0.7)
  ) +
  scale_shape_manual(
    name = "Direction",
    values = c("Higher in TR" = 16, "Higher in SW" = 17),
    guide = guide_legend(override.aes = list(size = 5)))+
  labs(
    title = "Top 20 Trait-SNPs with Largest Frequency Differences",
    x = "Associated Trait",
    y = "Variant ID"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

trait_compare

ggsave(trait_compare,
       filename = "/Users/xiyas/V2_Genome_reporting/Plots/trait_compare.pdf",width = 897/72,height = 530/72)


########## not comparing differencies, but rather the top carrier finidngs =====================
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(viridis)
library(patchwork)
library(tidytext)

# First, rename the populations for clarity and consistency
variant_data <- variant_freq %>%
  mutate(Population = ifelse(Population == "Turkish", "TR", 
                             ifelse(Population == "Swedish", "SW", Population)))

# For the first plot, get top variants for each population separately
# Get top 10 variants for TR
tr_top_vars <- variant_data %>%
  filter(Population == "TR") %>%
  arrange(desc(freq)) %>%
  head(10) %>%
  pull(variant_genotype)

# Get top 10 variants for SW
sw_top_vars <- variant_data %>%
  filter(Population == "SW") %>%
  arrange(desc(freq)) %>%
  head(10) %>%
  pull(variant_genotype)

# Combine the lists and get the unique variants
combined_top_vars <- unique(c(tr_top_vars, sw_top_vars))

# Filter for only those top variants
viz_data <- variant_data %>%
  filter(variant_genotype %in% combined_top_vars) %>%
  # Create a combined label for better readability
  mutate(var_display_name = paste0(variant_genotype, "\n", Traits.name))

# Create the faceted bar plot with TR and SW sorted separately
bar_plot <- viz_data %>%
  # This ensures each facet is sorted independently
  ggplot(aes(x = reorder_within(var_display_name, freq, Population), y = freq, fill = Population)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("SW" = "#3b9ab2", "TR" = "#e1812c")) +
  # The scale_x_reordered() function works with reorder_within()
  scale_x_reordered() +
  facet_wrap(~ Population, scales = "free_y") +
  labs(title = "Top Genetic Variants by Frequency",
       x = NULL,
       y = "Carrier Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  coord_flip()


bar_plot
# 2. Heatmap of variants by population
heat_plot <- viz_data %>%
  ggplot(aes(x = Population, y = fct_reorder(var_display_name, freq), fill = freq)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Frequency") +
  labs(title = "Variant Frequency Heatmap",
       x = "Population",
       y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")

# 3. Scatter plot comparing frequencies between populations
# First reshape the data to have TR and SW side by side
pop_comparison <- variant_data %>%
  select(variant_genotype, Traits.name, Population, freq, category) %>%
  pivot_wider(names_from = Population, values_from = freq, values_fill = 0) %>%
  filter(!is.na(TR) & !is.na(SW))  # Keep only variants present in both populations

scatter_plot <- ggplot(pop_comparison, aes(x = TR, y = SW, color = category)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_viridis_d(option = "turbo") +
  xlim(0, max(c(pop_comparison$TR, pop_comparison$SW)) * 1.05) +
  ylim(0, max(c(pop_comparison$TR, pop_comparison$SW)) * 1.05) +
  labs(title = "Variant Frequency: TR vs SW",
       x = "TR Population Frequency",
       y = "SW Population Frequency",
       color = "Category") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")
scatter_plot
# 4. Top variants by category for each population
bubble_plot <- variant_data %>%
  group_by(category, Population) %>%
  summarise(max_freq = max(freq),
            median_freq = median(freq),
            variant_count = n_distinct(variant_genotype),
            .groups = "drop") %>%
  ggplot(aes(x = fct_reorder(category, median_freq), y = median_freq, 
             fill = Population, size = variant_count)) +
  geom_point(shape = 21, alpha = 0.8) +
  scale_fill_manual(values = c("SW" = "#3b9ab2", "TR" = "#e1812c")) +
  scale_size_continuous(range = c(3, 12)) +
  labs(title = "Median Frequency by Category",
       x = NULL,
       y = "Median Carrier Frequency",
       size = "Number of\nVariants",
       fill = "Population") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  coord_flip()

# 5. Distribution plot by population
ridgeline_plot <- ggplot(variant_data, aes(x = freq, fill = Population)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("SW" = "#3b9ab2", "TR" = "#e1812c")) +
  facet_wrap(~category, scales = "free_y") +
  labs(title = "Frequency Distribution by Category and Population",
       x = "Variant Frequency",
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top")

# Combine plots using patchwork
multi_panel_plot <- (bar_plot | heat_plot) / (scatter_plot | bubble_plot)
multi_panel_plot + plot_annotation(
  title = "Genetic Variant Frequency Comparison Between TR and SW Populations",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
)

# For individual plots
# bar_plot     # Top variants bar chart
# heat_plot    # Heatmap
# scatter_plot # Population comparison
# bubble_plot  # Category bubble chart
# ridgeline_plot # Distribution by category