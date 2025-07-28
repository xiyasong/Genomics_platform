library(dplyr)
public_AF_data <- bind_rows(
  clinvar_TR %>% select(Variant_info, MAX_AF, Population,ClinVar_CLNSIG,Database,ReviewStar,MAX_AF_Category),
  clinvar_SW %>% select(Variant_info, MAX_AF, Population,ClinVar_CLNSIG,Database,ReviewStar,MAX_AF_Category)
) %>%
  mutate(public_AF = ifelse(is.na(MAX_AF), 0, MAX_AF)) %>%
  distinct(Variant_info, Population, .keep_all = TRUE)

# First get TOTAL patient counts for each population
total_patients_TR <- n_distinct(clinvar_TR$patientID)  # All Turkish patients
total_patients_SW <- n_distinct(clinvar_SW$patientID)  # All Swedish patients

# Calculate allele frequencies
variant_level_freq <- bind_rows(
  # Turkish cohort
  clinvar_TR %>%
    group_by(Variant_info) %>%
    summarise(
      allele_count = sum(ifelse(Zygosity == "homo", 2, 1)),  # Homo=2, Hetero=1
      inhouse_AF = allele_count / (2 * total_patients_TR),    # 2N denominator
      Population = "Turkish",
      .groups = "drop"
    ),
  
  # Swedish cohort
  clinvar_SW %>%
    group_by(Variant_info) %>%
    summarise(
      allele_count = sum(ifelse(Zygosity == "homo", 2, 1)),
      inhouse_AF = allele_count / (2 * total_patients_SW),
      Population = "Swedish",
      .groups = "drop"
    )
)

final_comparison <- variant_level_freq %>%
  left_join(
    public_AF_data %>% select(Variant_info, Population, public_AF,ClinVar_CLNSIG,Database,ReviewStar,MAX_AF_Category),
    by = c("Variant_info", "Population")
  ) %>%
  mutate(
    AF_ratio = inhouse_AF / public_AF,
    
    AF_discrepancy = case_when(
      # 处理public_AF=0的特殊情况
      public_AF == 0 ~ "Novel variant",
      
      # 极端差异（优先判断）
      AF_ratio > 3 ~ "In-house >3x public",
      AF_ratio < 0.33 ~ "In-house <0.33x public",
      
      # 中间梯度划分
      AF_ratio >= 1.2 ~ "In-house 1.2-3x public",
      AF_ratio <= 0.8 ~ "In-house 0.33-0.8x public",
      
      # 核心一致区间
      between(AF_ratio, 0.8, 1.2) ~ "Consistent (0.8-1.2x)",
      
      # 安全兜底（理论上不应触发）
      TRUE ~ "Uncategorized"
    )
  )


final_comparison<-  final_comparison%>%
  mutate(ClinVar_CLNSIG = case_when(
    grepl("Conflicting", ClinVar_CLNSIG, ignore.case = TRUE) ~ "CPLP",
    tolower(ClinVar_CLNSIG) == "pathogenic" ~ "P",
    tolower(ClinVar_CLNSIG) == "likely_pathogenic" ~ "LP",
    TRUE ~ ClinVar_CLNSIG
  ))
final_comparison <- final_comparison %>%
  mutate(Population = factor(Population, 
                             levels = c("Turkish", "Swedish"),  # Turkish first
                             ordered = TRUE))
# ggplot(final_comparison, aes(x = public_AF, y = inhouse_AF, color = ClinVar_CLNSIG)) +
#   geom_jitter(
#     aes(color = ClinVar_CLNSIG),
#     width = 0.05,   # Horizontal jitter (adjust based on your x-axis scale)
#     height = 0.05,  # Vertical jitter (adjust for y-axis)
#     alpha = 0.6,
#     size = 2
#   ) +
#   geom_abline(slope = 1, linetype = "dashed") +
#   scale_x_log10(labels = scales::scientific) + 
#   scale_y_log10(labels = scales::scientific) +
#   scale_color_jco() + 
#   facet_wrap(~Population) +
#   labs(
#     title = "In-house vs. Public Allele Frequencies",
#     x = "Public MAX_AF",
#     y = "In-house Allele Frequency",
#     color = "Pathogenicity"
#   )+theme_bw()
final_comparison_sorted <- final_comparison %>% 
  arrange(desc(inhouse_AF))
library(writexl)
# Save to Excel
write_xlsx(
  final_comparison_sorted,
  path = "/Users/xiyas/V2_Genome_reporting/data/final_comparison_sorted.xlsx"  # File path
)
#####ploting
Inhouse_vs_public <- ggplot(final_comparison, aes(x = public_AF, y = inhouse_AF)) +
  geom_jitter(
    aes(shape = ClinVar_CLNSIG, color = AF_discrepancy),  # 颜色映射AF差异，形状映射ClinVar
    width = 0.05,
    height = 0.05,
    alpha = 0.6,
    size = 3
  ) +
  geom_abline(slope = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10(labels = scales::scientific) + 
  scale_y_log10(labels = scales::scientific) +
  # 手动设置颜色映射（红蓝渐变）
  scale_color_manual(
    name = "AF Discrepancy",
    values = c(
      "In-house >3x public" = "#D73027",  # 深红（最高AF差异）
      "In-house 1.2-3x public" = "#FC8D59",  # 橙红
      "Consistent (0.8-1.2x)" = "#FEE090",  # 浅黄/中性色
      "In-house 0.33-0.8x public" = "#91BFDB",  # 浅蓝
      "In-house <0.33x public" = "#313695",  # 深蓝（最低AF差异）
      "Novel variant" = "#808080"  # 灰色（特殊类别）
      #[1] "#A50026" "#D73027" "#F46D43" "#FDAE61" "#FEE090" "#E0F3F8" "#ABD9E9" "#74ADD1" "#4575B4" "#313695"
    ),
    breaks = c(  # 控制图例顺序
      "In-house >3x public",
      "In-house 1.2-3x public",
      "Consistent (0.8-1.2x)",
      "In-house 0.33-0.8x public",
      "In-house <0.33x public",
      "Novel variant"
    )
  ) +
  # 形状映射（ClinVar分类）
  scale_shape_manual(
    name = "ClinVar Significance",
    values = c(
      "P" = 16,    # 实心圆
      "LP" = 17,   # 实心三角
      "CPLP" = 4   # 叉号
    ),
    labels = c(
      "P" = "Pathogenic",
      "LP" = "Likely Pathogenic",
      "CPLP" = "Conflicting"
    )
  ) +
  facet_grid(
      Population ~  ReviewStar ) +  # 按Population分面
  labs(
    title = "In-house vs. Public Allele Frequencies (facet by Review Star 0,1,2,3,4)",
    x = "Public MAX_AF (log10)",
    y = "In-house AF (log10)"
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
    
  )

# 显示图形
print(Inhouse_vs_public)
ggsave(Inhouse_vs_public,
       filename = "/Users/xiyas/V2_Genome_reporting/Plots/Inhouse_vs_public_neweset.pdf",width = 12,height =8)
       

# # Reviewstar 做法 ===================
# final_comparison$ReviewStar <- as.character(final_comparison$ReviewStar)
# Inhouse_vs_public<- ggplot(final_comparison, aes(x = public_AF, y = inhouse_AF)) +
#   geom_jitter(
#     aes(shape = ClinVar_CLNSIG,color = ReviewStar),
#     width = 0.05,   # Horizontal jitter (adjust based on your x-axis scale)
#     height = 0.05,  # Vertical jitter (adjust for y-axis)
#     alpha = 0.6,
#     size = 2
#   ) +
#   scale_color_nejm()+
#   geom_abline(slope = 1, linetype = "dashed") +
#   # How to Verify the Line is Correct
#   #annotate("point", x = 1e-05, y = 1e-05, color = "red", size = 3) +  # 1e-05 marker
#   #annotate("point", x = 1e-03, y = 1e-03, color = "red", size = 3) +  # 1e-03 marker
#   scale_x_log10(labels = scales::scientific) + 
#   scale_y_log10(labels = scales::scientific) +
#   scale_shape_manual(values = c(
#     "P" = 16,        # Solid circle
#     "LP" = 17, # Triangle
#     "CPLP" = 4       # X mark
#   )) +
#   facet_wrap(~Population) +
#   labs(
#     title = "In-house vs. Public Max_AF",
#     x = "Public MAX_AF",
#     y = "In-house Allele Frequency",
#     color = "ReviewStar"
#   )+theme_bw(base_size = 18)
# 
# ggsave(Inhouse_vs_public,
#        filename = "/Users/xiyas/V2_Genome_reporting/Plots/Inhouse_vs_public.pdf",width = 752/72,height = 443/72)

############ backup：老的6子图做法 ===============
# # 确保分类变量已转换为因子，并设定顺序
# final_comparison <- final_comparison %>%
#   mutate(
#     Population = factor(Population, levels = c("Turkish", "Swedish")),
#     ClinVar_CLNSIG = factor(ClinVar_CLNSIG, levels = c("P", "LP", "CPLP")),
#     AF_discrepancy = factor(AF_discrepancy),
#     ReviewStar = factor(ReviewStar, levels = c("0","1","2","3","4"))
#   )
# 
# # 定义形状映射（AF_discrepancy）
# shape_values <- c(
#   "Consistent" = 16,          # 实心圆
#   "In-house <0.5x public" = 17, # 实心三角
#   "In-house >2x public" = 15,   # 实心方形
#   "Novel variant" = 18          # 实心菱形
# )
# 
# # 分面组合：Population * ClinVar_CLNSIG
# Inhouse_vs_public_6_plot <- ggplot(final_comparison, aes(x = public_AF, y = inhouse_AF)) +
#   geom_jitter(
#     aes(shape = AF_discrepancy, color = ReviewStar),  # 形状=AF差异，颜色=ReviewStar
#     width = 0.05,
#     height = 0.05,
#     alpha = 0.7,
#     size = 3.5
#   ) +
#   geom_abline(
#     slope = 1, 
#     linetype = "dashed", 
#     color = "gray40", 
#     linewidth = 0.8
#   ) +
#   scale_x_log10(
#     labels = scales::scientific,
#   ) + 
#   scale_y_log10(
#     labels = scales::scientific,
#   ) +
#   scale_color_nejm(
#     name = "Review Star",  # 颜色图例标题
#     labels = c("0 (Low)", "1", "2", "3", "4 (High)")  # 自定义标签
#   ) +
#   scale_shape_manual(
#     name = "AF Discrepancy",  # 形状图例标题
#     values = shape_values,
#     guide = guide_legend(override.aes = list(size = 3))  # 调整图例符号大小
#   ) +
#   facet_grid(
#     Population ~ ClinVar_CLNSIG,  # 分面：行=Population，列=ClinVar_CLNSIG
#     labeller = labeller(  # 自定义分面标签
#       Population = c("Turkish" = "Turkish", "Swedish" = "Swedish"),
#       ClinVar_CLNSIG = c("P" = "Pathogenic", "LP" = "Likely Pathogenic", "CPLP" = "Conflicting")
#     )
#   ) +
#   labs(
#     title = "In-house vs. Public Allele Frequencies by Population and ClinVar Significance",
#     x = "Public MAX_AF",
#     y = "In-house Allele Frequency"
#   ) +
#   theme_bw(base_size = 22) +
#   theme(
#     legend.position = "right",
#     legend.box = "vertical",
#     legend.spacing.y = unit(0.5, "cm"),
#     strip.text = element_text(face = "bold", size = 12),  # 分面标题加粗
#     strip.background = element_rect(fill = "white"),      # 分面标题背景
#     panel.grid.minor = element_blank()  # 隐藏次要网格线
#   )
# Inhouse_vs_public_6_plot
# ggsave(Inhouse_vs_public_6_plot,
#        filename = "/Users/xiyas/V2_Genome_reporting/Plots/Inhouse_vs_public_6_plot.pdf",width = 1229/72,height = 675/72)
