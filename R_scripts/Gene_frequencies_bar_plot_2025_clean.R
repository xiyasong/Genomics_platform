# 定义 cohort 变量（假设 cohort 是用户选择的输入）
cohort <- "TR"  # 或者 "SW"

# 根据 cohort 选择数据表
if (cohort == "TR") {
  # TR 数据
  temp_table <- sorted_tab_TR_filter$sorted_tab_old_filter
  temp_table_sillico <- sorted_tab_TR_sillico$sorted_tab_old_sillico
} else if (cohort == "SW") {
  # SW 数据
  temp_table <- sorted_tab_SW_filter$sorted_tab_old_filter
  temp_table_sillico <- sorted_tab_SW_sillico$sorted_tab_old_sillico
} else {
  # 如果 cohort 不是 TR 或 SW，抛出错误或设置默认值
  stop("Invalid cohort value. Please choose 'TR' or 'SW'.")
}

## genes that needed for ClinVar figure =============
gene_counts_TR <- temp_table %>%
  group_by(Genes) %>%
  summarise(Total_Freq = sum(Freq)) %>% arrange(desc(Total_Freq))

gene_selected_TR = gene_counts_TR$Genes[1:30]

# genes that needed for sillico figure ==============
gene_counts_TR_silloco <- temp_table_sillico %>%
  group_by(Genes) %>%
  summarise(Total_Freq = sum(Freq)) %>% arrange(desc(Total_Freq))

gene_selected_TR_sillico = gene_counts_TR_silloco$Genes[1:30]


##  dataframe that needed for ClinVar figure =============
temp_table <- 
  temp_table %>% 
  mutate(group = case_when(Freq == 1 ~ "count(1)",
                           Freq > 1 & Freq <= 5 ~ "count(2~5)",
                           Freq > 5 & Freq <= 10 ~ "count(5~10)",
                           Freq > 10 & Freq <= 15 ~ "count(10~15)",
                           Freq > 15 & Freq <= 20 ~ "count(15~20)",
                           Freq > 20 & Freq <= 25 ~ "count(20~25)",
                           Freq > 25 & Freq <= 30 ~ "count(25~30)",
                           Freq > 30 & Freq <= 35 ~ "count(30~35)",
                           Freq > 35 & Freq <= 40 ~ "count(35~40)",
                           Freq > 40 ~ "count(>40)"
  )) %>% 
  arrange(-desc(Genes), -desc(group)) %>% 
  as_tibble()
# Selected df with top30 genes
temp_table_selected <- 
  temp_table[temp_table$Genes %in%gene_selected_TR,] %>% 
  arrange(-desc(Genes), desc(Freq))
desired_levels <- c(
  "count(1)", "count(2~5)", "count(5~10)", "count(10~15)", "count(15~20)",
  "count(20~25)", "count(25~30)", "count(30~35)", "count(35~40)", "count(>40)"
)

temp_table_selected$group <- factor(temp_table_selected$group,levels = desired_levels)

##  dataframe that needed for sillico figure =============
# adding frequency group information
temp_table_sillico <- 
  temp_table_sillico %>% 
  mutate(group = case_when(Freq == 1 ~ "count(1)",
                           Freq > 1 & Freq <= 5 ~ "count(2~5)",
                           Freq > 5 & Freq <= 10 ~ "count(5~10)",
                           Freq > 10 & Freq <= 15 ~ "count(10~15)",
                           Freq > 15 & Freq <= 20 ~ "count(15~20)",
                           Freq > 20 & Freq <= 25 ~ "count(20~25)",
                           Freq > 25 & Freq <= 30 ~ "count(25~30)",
                           Freq > 30 & Freq <= 35 ~ "count(30~35)",
                           Freq > 35 & Freq <= 40 ~ "count(35~40)",
                           Freq > 40 ~ "count(>40)"
  )) %>% 
  arrange(-desc(Genes), -desc(group)) %>% 
  as_tibble()

# Selected df with top30 ClinVar genes
temp_table_sillico_selected <- 
  temp_table_sillico[temp_table_sillico$Genes %in% gene_selected_TR,] %>% 
  arrange(-desc(Genes), desc(Freq))
temp_table_sillico_selected$Freq <- -(temp_table_sillico_selected$Freq)


#merge the two df to a merged_df
merge_col <- c('Genes','Freq','group','IMPACT')

#cut is just a cut of undeeded column values
temp_table_sillico_selected_cut <-temp_table_sillico_selected %>% select(merge_col)
temp_table_selected_cut <-temp_table_selected %>% select(merge_col)


# ploting way 2 for ClinVar Variants--- 
## combine the clinvar and sillico variants ================

temp_table_selected_1 <- 
  temp_table%>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>% 
  summarise(count = n()) %>% filter(Genes %in% gene_selected_TR) 

temp_table_sillico_selected_1 <- 
  temp_table_sillico %>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>% 
  summarise(count = n()) %>% filter(Genes %in% gene_selected_TR)


merge_df_TR <- left_join(temp_table_selected_1,temp_table_sillico_selected_1,by=c("Genes","group"))
merge_df_TR[is.na(merge_df_TR)] <-0

merge_df_TR_long <- melt(merge_df_TR,id.vars=c("Genes","group"))
merge_df_TR_long[which(merge_df_TR_long$variable =="count.y"),]$value <-
  -(merge_df_TR_long[which(merge_df_TR_long$variable =="count.y"),]$value)

# graph setting
confidence_setting <- ifelse(gene_selected_TR %in% High_genes$Genes, "dark blue", "black")
merge_df_TR_long$group <- factor(merge_df_TR_long$group,levels = desired_levels)
## label the putative variants where impact == "HIGH"
#merge_df_TR$hLoFs_setting <- ifelse(merge_df_TR$Freq < 0 & merge_df_TR$IMPACT == "HIGH", TRUE, FALSE)
#text_data <- subset(merge_df_TR,hLoFs_setting ==TRUE)
## ploting way 2 ggplot code
Figure5D <-merge_df_TR_long %>% 
  ggplot(aes(x = Genes, y = value, fill = group)) +
  geom_bar(stat="identity", color = "white",alpha=0.9) +
  ####sort the x axis by the order using scale_x_discrete()!!
  scale_x_discrete(limits = gene_selected_TR) +
  scale_y_continuous(expand = c(0,0), limits = c(-15, 15)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,8,10,5,4,3,1)]) +
  ggtitle("SW Carrier frequency of ClinVar P/LP variants") + 
  geom_hline(yintercept=0)+
  ylab("No.unique variant") + 
  geom_text(data=subset(merge_df_TR_long,value !=0),mapping = aes(label = abs(value)), position = position_stack(vjust = 0.5))+
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_line(size = 0.5, colour = "grey"),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_markdown(size = 12, colour = "black", angle = 315, vjust = 0.5, hjust = 0,color = confidence_setting)
  )

Figure5D
#ggsave(Figure5D,filename = "/Users/xiyas/V2_Genome_reporting/Plots/Figure5D.pdf",width = 9,height = 6)


## sillico plot ===================

temp_table_sillico_selected <- 
  temp_table_sillico %>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>% 
  summarise(count = n()) %>% filter(Genes %in% gene_selected_TR_sillico)

temp_table_selected <- 
  temp_table%>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>% 
  summarise(count = n()) %>% filter(Genes %in% gene_selected_TR_sillico) 

merge_df_TR <- left_join(temp_table_sillico_selected,temp_table_selected,by=c("Genes","group"))
merge_df_TR[is.na(merge_df_TR)] <-0

merge_df_TR_long <- melt(merge_df_TR,id.vars=c("Genes","group"))
merge_df_TR_long[which(merge_df_TR_long$variable =="count.y"),]$value <- -(merge_df_TR_long[which(merge_df_TR_long$variable =="count.y"),]$value)


# Newest version plot : TR -----------------------------------
High_genes <-GeneDB %>% filter(Gene.Disease.confidence.level == "High_confidence") %>% select(Genes)

confidence_setting <- ifelse(gene_selected_TR_sillico %in% High_genes$Genes, "dark blue", "black")
merge_df_TR_long$group <- factor(merge_df_TR_long$group,levels = desired_levels)
if (cohort == "TR") {
  fill_colors <- RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,7,8,9,10,5,4,3,2,1)]
  plot_title <- "TR Carrier frequency of pLOFs and p-risk variants"
} else if (cohort == "SW") {
  fill_colors <- RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,7,8,9,5,1)]
  plot_title <- "SW Carrier frequency of pLOFs and p-risk variants"
} else {
  stop("Invalid cohort value. Please choose 'TR' or 'SW'.")
}

# 绘制图表
sillico_plot <- merge_df_TR_long %>% 
  ggplot(aes(x = Genes, y = value, fill = group)) +
  geom_bar(stat = "identity", color = "white", alpha = 0.9) +
  scale_x_discrete(limits = gene_selected_TR_sillico) +
  scale_y_continuous(expand = c(0, 0), limits = c(-6, 50)) +
  theme_bw() +
  scale_fill_manual(values = fill_colors) +  # 动态设置颜色
  ggtitle(plot_title) +  # 动态设置标题
  geom_hline(yintercept = 0) +
  ylab("No. unique variant") +
  geom_text(
    data = subset(merge_df_TR_long, value != 0),
    mapping = aes(label = abs(value)), 
    position = position_stack(vjust = 0.5)
  ) +
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_line(size = 0.5, colour = "grey"),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_markdown(
      size = 12, 
      colour = "black", 
      angle = 315, 
      vjust = 0.5, 
      hjust = 0,
      color = confidence_setting  # 动态设置 x 轴文本颜色
    )
  )

# 显示图表
print(sillico_plot)
sillico_plot

#ggsave(sillico_plot,filename = "/Users/xiyas/V2_Genome_reporting/Plots/Figure5F.pdf",width = 9,height = 6)

