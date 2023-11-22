library(dplyr)
library(ggsci)
library(ggpubr)
library(ggtext)
library(tidyverse)
library(ggsci)
library(reshape2)

####new stacked barplot, based on the filtered results -----------
df_clinvar_filter <- df_clinvar %>% filter(MAX_AF<0.05)%>%filter(ClinVar_CLNSIG!='Conflicting_interpretations_of_pathogenicity')
df_puta_filter <- df_puta

#For the df_puta_filter --------------------------------
test<- as.data.frame(table(df_puta_filter$Variant_info))
colnames(test)[1]<- "Variant_info"
#df_clinvar_filter_Summary <- merge(test,df_clinvar_filter,by = "SZAID")
df_puta_filter_Summary <- merge(test,df_puta_filter,by = "Variant_info")
df_puta_filter_Summary <- df_puta_filter_Summary[!duplicated(df_puta_filter_Summary$Variant_info), ]
df_puta_filter_Summary <- 
  df_puta_filter_Summary %>% 
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

df_puta_filter_Summary_list <- 
  df_puta_filter_Summary %>% 
  dplyr::select(Freq, Genes) %>% 
  group_by(Genes) %>% 
  summarise(Freq = sum(Freq)) %>% 
  arrange(desc(Freq))

df_puta_filter_gene_selected <- 
  df_puta_filter_Summary_list$Genes[1:30]

df_puta_filter_Summary_selected <- 
  df_puta_filter_Summary[df_puta_filter_Summary$Genes %in% df_puta_filter_gene_selected,] %>% 
  arrange(-desc(Genes), desc(Freq))
 
df_puta_filter_Summary1 <- 
  df_puta_filter_Summary %>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>% 
  summarise(count = n()) %>% filter(Genes %in% df_puta_filter_gene_selected)

##add a new count from clinvar classes
new_count<- df_clinvar_filter_Summary %>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>%
  summarise(count = n()) %>% filter(Genes %in% df_puta_filter_gene_selected) 

df_puta_filter_Summary1 <- left_join(df_puta_filter_Summary1,new_count,by=c("Genes","group"))
df_puta_filter_Summary1[is.na(df_puta_filter_Summary1)] <-0

###reshape
df_puta_filter_Summary1_long <- melt(df_puta_filter_Summary1,id.vars=c("Genes","group"))
df_puta_filter_Summary1_long[which(df_puta_filter_Summary1_long$variable =="count.y"),]$value <- -(df_puta_filter_Summary1_long[which(df_puta_filter_Summary1_long$variable =="count.y"),]$value)

High_genes <-GeneDB %>% filter(Gene.Disease.confidence.level == "High_confidence") %>% select(Genes)
confidence_setting <- ifelse(df_puta_filter_gene_selected %in% High_genes$Genes, "dark blue", "black")
###for df_puta:
all_levels <- c(
  "count(1)", "count(2~5)", "count(5~10)", "count(10~15)", "count(15~20)",
  "count(20~25)", "count(25~30)", "count(30~35)", "count(35~40)", "count(>40)"
)
# Get unique levels from the selected data frame

df_puta_filter_Summary_selected$group <-
  factor(df_puta_filter_Summary_selected$group, levels =all_levels,ordered=TRUE)
df_puta_filter_Summary1$group <-
  factor(df_puta_filter_Summary1$group, levels = all_levels,ordered = TRUE)
df_puta_filter_Summary1_long$group <-
  factor(df_puta_filter_Summary1_long$group, levels = all_levels,ordered=TRUE)
# Newest version plot -----------------------------------
df_puta_filter_Summary1_long %>% 
  ggplot(aes(x = Genes, y = value, fill = group)) +
  geom_bar(stat="identity", color = "white",alpha=0.9) +
  ####sort the x axis by the order using scale_x_discrete()!!
  scale_x_discrete(limits = df_puta_filter_gene_selected) +
  scale_y_continuous(expand = c(0,0), limits = c(-6, 50)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,7,8,9,10,5,4,3,1)]) +
  ggtitle("Turkish putative carrier frequency") + 
  geom_hline(yintercept=0)+
  geom_text(data=subset(df_puta_filter_Summary1_long,value !=0),mapping = aes(label = abs(value)), position = position_stack(vjust = 0.5))+
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
#ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot.pdf",height = 5,width =9)
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_puta_turkish_207.pdf",height = 5,width =9)

# {not used}ggplot2 for stacked plot which with variant total numbers-----------------------------------------------------------------

df_puta_filter_Summary_selected %>%  
  ggplot(aes(x = Genes, y = Freq, fill = group)) +
  geom_bar(stat="identity", color = "white") +
  scale_x_discrete(limits = df_puta_filter_gene_selected) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 350)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,7,8,9,10,5,4,3,2,1)]) +
  ggtitle("test") +
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_line(size = 0.5, colour = "grey"),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 315, vjust = 0.5, hjust = 0)
  )

##add a new count from puta
new_count<- df_puta_filter_Summary %>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>%
  summarise(count = n()) %>% filter(Genes %in% df_puta_filter_gene_selected) 

df_puta_filter_Summary1 <- left_join(df_puta_filter_Summary1,new_count,by=c("Genes","group"))
df_puta_filter_Summary1[is.na(df_puta_filter_Summary1)] <-0


###reshape
df_clinvar_filter_Summary1_long <- melt(df_clinvar_filter_Summary1,id.vars=c("Genes","group"))
df_clinvar_filter_Summary1_long[which(df_clinvar_filter_Summary1_long$variable =="count.y"),]$value <- -(df_clinvar_filter_Summary1_long[which(df_clinvar_filter_Summary1_long$variable =="count.y"),]$value)

###adding confidence info
High_genes <-GeneDB %>% filter(Gene.Disease.confidence.level == "High_confidence") %>% select(Genes)
confidence_setting <- ifelse(gene_selected %in% High_genes$Genes, "dark blue", "black")

# V1:ggplot2 for stacked plot which with unique variants-----------------------------------
df_puta_filter_Summary1_long %>% 
  ggplot(aes(x = Genes, y = value, fill = group)) +
  geom_bar(stat="identity", color = "white",alpha=0.9) +
  scale_x_discrete(limits = df_puta_filter_gene_selected) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,7,8,9,10,5,4,3,2,1)]) +
  ggtitle("test") + 
  geom_text(mapping = aes(label = value), position = position_stack(vjust = 0.5))+
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
#ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot.pdf",height = 5,width =9)
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot_df_puta_filter.pdf",height = 5,width =9)

#V2:
High_genes <-GeneDB %>% filter(Gene.Disease.confidence.level == "High_confidence") %>% select(Genes)
confidence_setting <- ifelse(df_puta_filter_gene_selected %in% High_genes$Genes, "dark blue", "black")

df_puta_filter_Summary1_long %>% 
  ggplot(aes(x = Genes, y = value, fill = group)) +
  geom_bar(stat="identity", color = "white",alpha=0.9) +
  scale_x_discrete(limits = df_puta_filter_gene_selected) +
  scale_y_continuous(expand = c(0,0), limits = c(-10, 35)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,7,8,9,10,5,4,3,2,1)]) +
  ggtitle("test") + 
  geom_hline(yintercept=0)+
  geom_text(data=subset(df_puta_filter_Summary1_long,value !=0),mapping = aes(label = abs(value)), position = position_stack(vjust = 0.5))+
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
#ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot.pdf",height = 5,width =9)
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot_df_puta_filter_v2.pdf",height = 5,width =9)


###For the df_clinvar_filter --------------
## you table with SZAID and table with Variant_info will show different results.

### sorted_tab_TR_filter$sorted_tab_old_filter = df_clinvar_filter_Summary
test<- as.data.frame(table(df_clinvar_filter$SZAID))

colnames(test)[1]<- "SZAID"
#df_clinvar_filter_Summary <- merge(test,df_clinvar_filter,by = "SZAID")
df_clinvar_filter_Summary <- merge(test,df_clinvar_filter,by = "SZAID")
df_clinvar_filter_Summary <- df_clinvar_filter_Summary[!duplicated(df_clinvar_filter_Summary$SZAID), ]

df_clinvar_filter_Summary <- 
  df_clinvar_filter_Summary %>% 
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

## top 30 genes ------
df_clinvar_filter_Summary_list <- 
  df_clinvar_filter_Summary %>% 
  dplyr::select(Freq, Genes) %>% 
  group_by(Genes) %>% 
  summarise(Freq = sum(Freq)) %>% 
  arrange(desc(Freq))

df_clinvar_filter_gene_selected <- 
  df_clinvar_filter_Summary_list$Genes[1:30]

df_clinvar_filter_Summary_selected <- 
  df_clinvar_filter_Summary[df_clinvar_filter_Summary$Genes %in% df_clinvar_filter_gene_selected,] %>% 
  arrange(-desc(Genes), desc(Freq))


###########After selected genes  -----
df_clinvar_filter_Summary1 <- 
  df_clinvar_filter_Summary %>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>% 
  summarise(count = n()) %>% filter(Genes %in% df_clinvar_filter_gene_selected)

##add a new count from putative variants
new_count<- df_puta_filter_Summary %>% 
  dplyr::select(Genes, group) %>% 
  group_by(Genes, group) %>%
  summarise(count = n()) %>% filter(Genes %in% df_clinvar_filter_gene_selected) 

df_clinvar_filter_Summary1 <- left_join(df_clinvar_filter_Summary1,new_count,by=c("Genes","group"))
df_clinvar_filter_Summary1[is.na(df_clinvar_filter_Summary1)] <-0

###for df_clinvar_filter:
all_levels <- c(
  "count(1)", "count(2~5)", "count(5~10)", "count(10~15)", "count(15~20)",
  "count(20~25)", "count(25~30)", "count(30~35)", "count(35~40)", "count(>40)"
)
# Get unique levels from the selected data frame

df_clinvar_filter_Summary_selected$group <-
  factor(df_clinvar_filter_Summary_selected$group, levels =all_levels,ordered=TRUE )

df_clinvar_filter_Summary1$group <-
factor(df_clinvar_filter_Summary1$group, levels = all_levels,ordered = TRUE)
###reshape
df_clinvar_filter_Summary1_long <- melt(df_clinvar_filter_Summary1,id.vars=c("Genes","group"))
df_clinvar_filter_Summary1_long[which(df_clinvar_filter_Summary1_long$variable =="count.y"),]$value <- -(df_clinvar_filter_Summary1_long[which(df_clinvar_filter_Summary1_long$variable =="count.y"),]$value)

###check  unique(df_clinvar_filter_Summary1_long$group)
df_clinvar_filter_Summary1_long$group <-
  factor(df_clinvar_filter_Summary1_long$group, levels = all_levels,ordered=TRUE)

# {not used}ggplot2 for stacked plot which with variant total numbers-----------------------------------------------------------------

df_clinvar_filter_Summary_selected %>%  
  ggplot(aes(x = Genes, y = Freq, fill = group)) +
  geom_bar(stat="identity", color = "white") +
  scale_x_discrete(limits = df_clinvar_filter_gene_selected) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 30)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,8,10,5,4,3,1)]) +
  ggtitle("Turkish P/LP carrier frequency") +
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_line(size = 0.5, colour = "grey"),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 315, vjust = 0.5, hjust = 0)
  )
#ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant2.pdf",height = 5,width =9)
new_select <- df_clinvar_filter_Summary1_long$Genes

#df_clinvar_filter_Summary1_long  <- df_clinvar_filter_Summary1_long[match(df_clinvar_filter_gene_selected, df_clinvar_filter_Summary1_long$Genes), ]  
# Newest version plot -----------------------------------
#ggplot2 for stacked plot which with unique variants

#V2:### The plot that adding confidence info for the x axis "blue" "dark" -------------
High_genes <-GeneDB %>% filter(Gene.Disease.confidence.level == "High_confidence") %>% select(Genes)
confidence_setting <- ifelse(df_clinvar_filter_gene_selected %in% High_genes$Genes, "dark blue", "black")

df_clinvar_filter_Summary1_long %>% 
  ggplot(aes(x = Genes, y = value, fill = group)) +
  geom_bar(stat="identity", color = "white",alpha=0.9) +
  ####sort the x axis by the order using scale_x_discrete()!!
  scale_x_discrete(limits = df_clinvar_filter_gene_selected) +
  scale_y_continuous(expand = c(0,0), limits = c(-10, 15)) +
  theme_bw() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")[c(6,8,10,5,4,3,1)]) +
  ggtitle("Turkish P/LP carrier frequency") + 
  geom_hline(yintercept=0)+
  geom_text(data=subset(df_clinvar_filter_Summary1_long,value !=0),mapping = aes(label = abs(value)), position = position_stack(vjust = 0.5))+
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
#ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot.pdf",height = 5,width =9)
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_clinvar_turkish_207.pdf",height = 5,width =9)

#ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot.pdf",height = 5,width =9)
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant_count_stacked_plot_df_puta_filter_v2.pdf",height = 5,width =9)

###secondary finding list

#########Frozen part -------------
df_clinvar_filter_gene <- as.data.frame(sort(table(Total_Clin_var_genes), decreasing=TRUE))
df_puta_filter_gene <- as.data.frame(sort(table(Total_puta_var_genes), decreasing=TRUE))

ggplot(df_clinvar_filter_gene, aes(Total_Clin_var_genes,Freq)) +geom_linerange(aes(x = (Total_Clin_var_genes, ymin = 0, ymax = Freq),
                                                                            color = "lightgray", size = 1.5)+
  geom_point(aes(color = Total_Clin_var_genes), size = 2)
  
  
  
##NEED INSTALL THE ggpubr package  
ggbarplot(df_clinvar_filter_gene%>%top_n(10), x = "Total_Clin_var_genes", y = "Freq",
          fill = "Total_Clin_var_genes",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = TRUE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_pubclean())+font("x.text", size = 8, vjust = 0.5)

ggdotchart(df_clinvar_filter_gene%>%top_n(10), x = "Total_Clin_var_genes", y = "Freq",
           color = "Total_Clin_var_genes",                                # Color by groups
           palette = "jco", # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "cyl",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(dfm$mpg),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )



new_df <-df_clinvar_filter_gene%>%top_n(30)
Confidence_total <- c()
for (row in 1:dim(new_df)[1])
{
  print(row)
  Confidence <- GeneDB[match(new_df[row,]$Total_Clin_var_genes,GeneDB$Genes),]$Gene.Disease.confidence.level
  Confidence_total <- append(Confidence_total,Confidence)
}
new_df$Confidence <- Confidence_total


ggdotchart(new_df, x = "Total_Clin_var_genes", y = "Freq",
           color = "Confidence",                                # Color by groups
           palette = "jco", # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           #group = "cyl",                                # Order by groups
           dot.size = 10,                                 # Large dot size
           label = round(new_df$Freq),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr(),)+ theme(axis.text.x = element_text(size = 15),
                                           axis.text.y = element_text(size = 15),
                                           axis.title.x = element_text(size = 20),
                                           axis.title.y = element_text(size = 20))
# ggplot2 theme
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_carrier_frequencies_clinvar.pdf",height = 13,width = 10)


new_df <-df_puta_filter_gene%>%top_n(30)
Confidence_total <- c()
for (row in 1:dim(new_df)[1])
{
  print(row)
  Confidence <- GeneDB[match(new_df[row,]$Total_puta_var_genes,GeneDB$Genes),]$Gene.Disease.confidence.level
  Confidence_total <- append(Confidence_total,Confidence)
}
new_df$Confidence <- Confidence_total
ggdotchart(new_df, x = "Total_puta_var_genes", y = "Freq",
           color = "Confidence",                                # Color by groups
           palette = "jco", # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           #group = "cyl",                                # Order by groups
           dot.size = 10,                                 # Large dot size
           label = round(new_df$Freq),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()+theme(axis.text.x = element_text(size = 15),
                                       axis.text.y = element_text(size = 15),
                                       axis.title.x = element_text(size = 20),
                                       axis.title.y = element_text(size = 20))                # ggplot2 theme
)
ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_carrier_frequencies_putative.pdf",height = 13,width = 10)
#########normal ggplot
ggplot(df_clinvar_filter_gene%>%top_frac(.1), aes(x=Total_Clin_var_genes, y=Freq)) + 
  geom_bar(stat = "identity") + theme_solarized() +  theme(legend.position = "none",axis.text.x = element_text(angle=45,vjust = 0.5,hjust = 0.5))



                                                    