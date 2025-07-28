##group pie plots based on unique patients ------------

#####Separate the AD and AR genes and count by Positive/Negative/Carrier ---------
groups <- c('Health predipositions','Carrier-screening','Newborn-screening','Heriditary-cancer risk syndrome')
count_total <- c()
KEY= 'Health predipositions'
KEY = 'Carrier-screening'
KEY= 'Health predipositions'
KEY= 'Heriditary-cancer risk syndrome'
KEY = "Newborn-screening"

#----------df_clinvar for screening results --------
for(KEY in groups){
  aa <- 
    filter(X.CHROM != 'chrM',
           MAX_AF < 0.05 | is.na(MAX_AF),
           ClinVar_CLNSIG %in% c('Pathogenic', 'Likely_pathogenic'),
           ReviewStar >= 2,
           acmg_classification %in% c('Pathogenic', 'Likely_pathogenic'),
           !(am_class %in% c('benign', 'likely_benign')))%>% filter(grepl(KEY,Target.group)) %>%
    mutate(
      condition = case_when(
        grepl("recessive", Inheritances, ignore.case = TRUE) & Zygosity == "Homozygous" ~ "Positive",
        grepl("recessive", Inheritances, ignore.case = TRUE) & Zygosity == "Heterozygous" ~ "Carrier",
        grepl("dominant", Inheritances, ignore.case = TRUE) & Zygosity == "Heterozygous" ~ "Positive",
        TRUE ~ "Unsure")
        
  aa_positive <- aa %>% filter(condition=="Positive")
  ############### one patient the results only count one time 
  count1 <- length(unique(aa%>% 
                            filter(grepl("Positive",condition)) %>%.$patientID))
  count2 <- length(unique(aa%>% 
                            filter(grepl("Carrier",condition)) %>%.$patientID))
  
  ####count3 = 总人数-positive-carrier
  count3 <- 101-count1-count2
  pp <- c(count1,count2,count3)
  count_total <- append(count_total,pp)
  #pie(c(count,275-count),labels = c("Carrier(n=39)","Non-Carrier(n=61)"),col = c("#F6AE2D", "#33658A"))
}

####value <- count_total 
test = data.frame(group <- c(rep('SF v3.1',3),rep('Carrier-screening',3),rep('Newborn-screening',3),
                               rep('Heriditary-cancer risk syndrome',3)
                               ),
                condition <- rep(c('Positive','Carrier','Negative'),4),
                value <- count_total)

test$subject <- factor(test$group)
test$credit <- factor(test$condition) 

mypal = pal_nejm("default", alpha = 0.9)(3)
mypal
mypal= c("#BC3C29E5","#E18727E5","#0072B5E5")
test$condition <- factor(test$condition,levels =c("Positive","Carrier","Negative"))

########figure of grouped pie chart --------
ggplot(data=test, aes(x=" ", y=value, group=condition, colour=condition, fill=condition)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ subject) +theme_void() +
  scale_fill_manual(values= mypal)+theme_bw()+
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.position = 'bottom')+
  ggtitle("Patient screening status based on applying different gene panels")

ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Groupped_bar_chart_turkish_207.pdf",height = 6,width =11)


######单个条形图 by gene panels #####
##########test plot of bar plot for patients in separate panels--------------
KEY= 'Health predipositions'
KEY = 'Carrier-screening'
KEY= 'Health predipositions'
KEY= 'Heriditary-cancer risk syndrome'
KEY = "Newborn-screening"

#target_genes <-unique(GeneDB %>% 
#                        filter(grepl(KEY,Target.group)) %>%.$Genes)

AD_genes <- 
  unique(Summary_DB %>% filter(grepl("Autosomal dominant",inheritances))%>%.$Genes)
AR_genes <- 
  unique(Summary_DB %>% filter(grepl("Autosomal recessive",inheritances))%>%.$Genes)
###start count the Positive Negative and Carrier -----------
###should not use genes, should use groups

aa_2 <- 
  df_clinvar_filter%>% filter(grepl(KEY,final_target_group)) %>%
  mutate(condition = case_when(Zygosity == "Homozygous" ~ "Positive",
                               Genes %in% AD_genes & Zygosity == "Heterozygous" ~ "Positive",
                               Genes %in% AR_genes & Zygosity == "Heterozygous" ~ "Carrier",
                               TRUE ~ "Negative"))
#aa_positive_2 <- aa_2 %>% filter(condition=="Positive")

aa_2 <-aa_2 %>% filter(!(condition == "Negative"))
aa_2$value <- rep(1)
aa_2 <- aa_2[!duplicated(aa[c("Genes","patientID")]),]

aa_2$paste.name <- paste(aa_2$Genes,"-",aa_2$Disease)



aa_2 %>% 
  ggplot(aes(x = reorder(paste.name,-value), y = value, fill = condition)) +
  geom_bar(stat="identity", color = "white",alpha=0.9) +
  #scale_x_discrete(limits = new_select) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = c("#FFB90F","#A52A2A")) +
  ggtitle("test") + 
  geom_hline(yintercept=0)+
  #geom_text(data=subset(aa,value !=0),mapping = aes(label = abs(value)), position = position_stack(vjust = 0.5))+
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_line(size = 0.5, colour = "grey"),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust =1,vjust = 0.5)
  )
aa_2

ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant.1.2.disease.pdf",height = 12,width =9)

######完全条形图 not limited by gene panels #####
##########test plot of bar plot for patients in separate panels--------------

disease_tab <- table(aa_2$Disease,aa_2$condition) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
colnames(sorted_tab) <- c('Disease','Freq')
#target_genes <-unique(GeneDB %>% 
#                        filter(grepl(KEY,Target.group)) %>%.$Genes)

aa_2 <- aa_2%>% filter(Disease %in% disease_tab$Var1[1:10]) %>% 
  mutate(condition = case_when(Zygosity == "Homozygous" ~ "Positive",
                               Genes %in% AD_genes & Zygosity == "Heterozygous" ~ "Positive",
                               Genes %in% AR_genes & Zygosity == "Heterozygous" ~ "Carrier",
                               TRUE ~ "Negative"))


aa_2 <-aa_2 %>% filter(!(condition == "Negative"))
aa_2$value <- rep(1)
aa_2 <- aa_2[!duplicated(aa[c("Genes","patientID")]),]

aa_2$paste.name <- paste(aa_2$Genes,"-",aa_2$Disease)


aa_2 %>% 
  ggplot(aes(x = reorder(paste.name,-value), y = value, fill = condition)) +
  geom_bar(stat="identity", color = "white",alpha=0.9) +
  #scale_x_discrete(limits = new_select) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  theme_bw() +
  # scale_fill_manual(values = c("#91D1C2", "#00A087", "#8491B4",
  #                              "#4DBBD5", "#E64B35")) +
  scale_fill_manual(values = c("#FFB90F","#A52A2A")) +
  ggtitle("test") + 
  geom_hline(yintercept=0)+
  #geom_text(data=subset(aa,value !=0),mapping = aes(label = abs(value)), position = position_stack(vjust = 0.5))+
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_line(size = 0.5, colour = "grey"),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust =1,vjust = 0.5)
  )


ggsave("/Users/xiyas/V2_Genome_reporting/Plots/Gene_variant.1.2.disease.pdf",height = 12,width =9)




plot_data <- df_temp_turkish %>%
  group_by(Target.group, condition) %>%
  summarize(count = n(), .groups = "drop") %>%
  # Make sure status has the desired factor order
  mutate(status = factor(condition, levels = c("Positive", "Carrier", "Unsure")))

# Define a custom color palette
mypal <- c("#BC3C29E5", "#E18727E5", "#0072B5E5")  # Red, Orange, Blue

# Create the grouped pie chart
ggplot(data = plot_data, aes(x = "", y = count, group = condition, fill = condition)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) + 
  facet_grid(. ~ Target.group) +
  scale_fill_manual(values = mypal) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.position = 'bottom',
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Patient screening status based on applying different gene panels",
    fill = "Status"
  )




# Function to calculate and prepare panel data
prepare_panel_data <- function(data, panel_name, total_patients = 275) {
  # Filter data for specific panel
  panel_data <- data %>% filter(grepl(panel_name, Target.group))
  
  # Print count of records for debugging
  cat("Panel:", panel_name, "- Total records:", nrow(panel_data), "\n")
  
  # Count unique patients with each status
  count_positive <- length(unique(panel_data %>% 
                                    filter(condition == "Positive") %>% 
                                    pull(patientID)))
  
  count_carrier <- length(unique(panel_data %>% 
                                   filter(condition == "Carrier") %>% 
                                   pull(patientID)))
  
  count_unsure <- length(unique(panel_data %>% 
                                  filter(condition == "Unsure") %>% 
                                  pull(patientID)))
  
  # Calculate negative as total minus positive minus carrier minus unsure
  count_negative <- total_patients - count_positive - count_carrier - count_unsure
  
  # Print all counts for debugging
  cat("  Positive:", count_positive, 
      "| Carrier:", count_carrier, 
      "| Negative:", count_negative, 
      "| Unsure:", count_unsure, "\n")
  
  # Return the counts
  return(c(count_positive, count_carrier, count_negative, count_unsure))
}

# Define your panel names
panels <- c("Health predipositions/Disease risk", "Carrier-screening", "Newborn-screening", "Heriditary-cancer risk syndrome")

# Initialize empty vector for counts
count_total <- c()

# Filter the dataframe once before the loop
df_temp_turkish_filtered_P_LP <- df_temp_turkish %>% 
  filter(X.CHROM != 'chrM',
         MAX_AF < 0.05 | is.na(MAX_AF),
         ClinVar_CLNSIG %in% c('Pathogenic', 'Likely_pathogenic'),
         ReviewStar >= 2,
         acmg_classification %in% c('Pathogenic', 'Likely_pathogenic'),
         !(am_class %in% c('benign', 'likely_benign')))


cat("SF v3.1 records after filtering:", nrow(sf_data_check), "\n")
# Process each panel
for(panel in panels) {
  counts <- prepare_panel_data(df_temp_turkish_filtered_P_LP, panel)
  count_total <- c(count_total, counts)
}

# Create the data frame for plotting - FIXED the comma missing between 'Negative' and 'Unsure'
test <- data.frame(
  group = rep(panels, each = 4),  # Changed from 3 to 4 for the four conditions
  condition = rep(c('Positive', 'Carrier', 'Negative', 'Unsure'), length(panels)),
  value = count_total
)

# Set up factor levels
test$subject <- factor(test$group)
test$condition <- factor(test$condition, levels = c("Positive", "Carrier", "Negative", "Unsure"))

# Define color palette
mypal <- c("#BC3C29E5", "#E18727E5", "#0072B5E5", "grey")  # Red, Orange, Blue, Grey

# Create the grouped pie chart
ggplot(data = test, aes(x = "", y = value, group = condition, fill = condition)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) + 
  facet_grid(. ~ subject) +
  theme_void() +
  scale_fill_manual(values = mypal) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.position = 'bottom',
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  ggtitle("Patient screening status based on applying different gene panels")

