########Reading python_output_files--Process to generate df_clinvar and df_puta ----------

###now it includes 275 files 
setwd("/Users/xiyas/V2_Genome_reporting/python_output_turkish_275/nodup_4_file/")
#generated nodup file 3 with score 12,13,14,15 and 7,8,9,10 variants-disease 
#setwd("/Users/xiyas/V2_Genome_reporting/python_output_file_swedish/")
files = list.files(path = "/Users/xiyas/V2_Genome_reporting/python_output_turkish_275/nodup_4_file/", pattern = "_4_nodup.txt")
#Mitochondrial genes removed
#temp = read.delim("/Users/xiyas/V2_Genome_reporting/python_output_file/P001_106.hard-filtered.vcf.gz_vep_annotated.vcf_outfile_sp_Inheritance_3_nodup.txt")list_clinvar = list()
list_clinvar = list()
list_puta = list()

for(i in files){
  temp = read.delim(file = i)
  # columns to paste together
  cols <- c( "X.CHROM", "POS" , "Genotype" )
  temp$POS <- gsub(" ","",temp$POS)
  # create a new column `Variant_info` with the three columns collapsed together
  temp$Variant_info <- apply(temp[,cols] , 1 , paste , collapse = "-" )
  Clin_var <-  temp %>% filter(grepl("ClinP_LP_var",Database_type))
  ###Only keep the unique variant for one samples ----------------------------This is the probelm, only show the first gene-disease matcing for this variant 
  ##New: this seems cause trouble when counting by gene panels
  Clin_var <- Clin_var[!duplicated(Clin_var$SZAID), ]
  ##Extract needed cols
  ##This time keep all the cols for checking
  Clin_var$patientID <- rep(i)
  colnames(Clin_var)[26] <- "Sample"
  ###combine the one person with another person: new method: dont' do rbind in for loop, rather create list
  list_clinvar[[i]] <- Clin_var
  puta_var <-  temp %>% filter(grepl("NovelTrans",SZAID))
  ###Only keep the unique variant for one samples
  puta_var <- puta_var[!duplicated(puta_var$SZAID),]
  #puta_var <- puta_var[,c('Genes','Variant_info',"Zygosity","ClinVar_CLNSIG","Consequence","MAX_AF","Gene.Disease.confidence.level","Target.group","Disease")]
  puta_var$patientID <- rep(i)
  colnames(puta_var)[26] <- "Sample"
  list_puta[[i]] <- puta_var
}

df_clinvar = do.call(rbind, list_clinvar)
df_clinvar$ClinVar_CLNSIG <- sapply(strsplit(df_clinvar$ClinVar_CLNSIG,"&"), `[`, 1)
df_clinvar$ClinVar_CLNSIG <- sapply(strsplit(df_clinvar$ClinVar_CLNSIG,"/"), `[`, 1)
df_clinvar <- df_clinvar %>% filter(X.CHROM != 'chrM')
df_puta = do.call(rbind, list_puta)
df_puta <- df_puta %>% filter(X.CHROM != 'chrM')


###group the two cohort
df_clinvar$Population <- "Turkish"
swe_df_clinvar$Population <- 'Swedish'
# Combine the two data frames into one
combined_df <- rbind(df_clinvar,swe_df_clinvar)
combined_df$Population <- factor(combined_df$Population, levels = c("Turkish", "Swedish"))

####per sample plot ====
patient_counts <- combined_df %>%
  group_by(ClinVar_CLNSIG, patientID,Population) %>%
  summarize(count = n()) %>%
  ungroup()

patient_counts <- patient_counts %>%
  mutate(ClinVar_CLNSIG = 
           ifelse(ClinVar_CLNSIG == "Conflicting_interpretations_of_pathogenicity", "Conflicts", ClinVar_CLNSIG))
summary_stats <- patient_counts %>%
  group_by(Population,ClinVar_CLNSIG) %>%
  summarize(
    min_count = min(count),
    max_count = max(count),
    median_count = median(count),
    mean_count = mean(count)
  )

### ====violin plot based on per sample statistics
# Create the violin plot
ggplot(patient_counts, aes(x = ClinVar_CLNSIG, y = count, fill = Population)) +
  geom_violin(scale='width',width=0.8) +
  #geom_dotplot(aes(x = ClinVar_CLNSIG, y = count), binaxis = "y", stackdir = "center", dotsize = 0.5, fill = "black") +  # Add dots
  geom_boxplot(width = 0.2, fill = "white", color = "black") +  # Add boxplots for clarity
  labs(x = "Pathogenicity", y = "Variant Counts per sample") +
  scale_fill_nejm()  +  # Customize fill colors
  theme_minimal()+
  facet_grid(.~Population, scales = "free_x", space = "free_x") +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # Adjust font sizes
    axis.text.x = element_text(size = 14),     # Font size for x-axis labels
    axis.title.x = element_text(size = 16),    # Font size for x-axis title
    axis.text.y = element_text(size = 14),     # Font size for y-axis labels
    axis.title.y = element_text(size = 16),    # Font size for y-axis title
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  geom_text(data = summary_stats, aes(x = ClinVar_CLNSIG, y= max_count,label = paste("Min:", min_count, "\nMax:", max_count, "\nMedian:", median_count)),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 3, fontface = "bold", color = "black", inherit.aes = FALSE)
  

####2: count the clinvar statictics on unique variants, get sorted tab(distinguish zygosity)=============
###get the no duplicated variants 
test <- df_clinvar[!duplicated(df_clinvar$SZAID),]

# star frome here=======================================================================================================================
# df_clinvar = clinvar_TR
# test = clinvar_TR_unique
# df_puta =sillico_pLoFs_TR

sorted_tab <- df_clinvar %>% group_by(SZAID, Zygosity)%>% summarise(count = n())%>%arrange(desc(count))
sorted_tab <- merge(x=sorted_tab,y=test,by='SZAID',all.x = TRUE)
sorted_tab <- sorted_tab %>% filter(X.CHROM != 'chrM')
sorted_tab <- sorted_tab %>% arrange(desc(count))
sorted_tab$ClinVar_CLNSIG <- sapply(strsplit(sorted_tab$ClinVar_CLNSIG,"&"), `[`, 1)

###get sorted tab_old (not distinguish zygosity)=============
sorted_tab_old <- table(df_clinvar$SZAID) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
colnames(sorted_tab_old) <- c('SZAID','Freq')
sorted_tab_old <- merge(x=sorted_tab_old,y=test,by='SZAID',all.x = TRUE)
sorted_tab_old <- sorted_tab_old %>% filter(X.CHROM != 'chrM')

###make a plot for variant summary count Figure 4A------------
noAF<- sorted_tab_old %>% filter(is.na(MAX_AF))
table(noAF$ClinVar_CLNSIG)
rareAF <- sorted_tab_old %>% filter(MAX_AF<0.05)
table(rareAF$ClinVar_CLNSIG)
CommonAF <- sorted_tab_old %>% filter(MAX_AF>=0.05)
table(CommonAF$ClinVar_CLNSIG)
#length(intersect(unique(sorted_tab_old$SZAID),unique(swe_sorted_tab_old$SZAID)))

######Stacked barchart ============

# create a dataset
Pathogenicity <- c(rep("Pathogenic" , 3) , rep("Likely Pathogenic" , 3) , rep("Conflicts" , 3))
Max_AF <- rep(c("No_public AF" , "Public AF < 0.05" , "Public AF > =0.05") , 3)
Variants_count <- c(40,274,21,18,66,1,2,342,64)
group <- c(rep("Turkish", 9))
data <- data.frame(Pathogenicity,Max_AF,Variants_count,group)
# combine with swedish--------
data <- rbind(data,data_swe)
data$group <- factor(data$group, levels = c("Turkish", "Swedish"))
group_sums <- aggregate(Variants_count ~ group, data, sum)
# Merge the sum of counts back into the data
data <- merge(data, group_sums, by = "group", suffixes = c("", "_sum"))
# Calculate the proportion within each group
data$Proportion <- data$Variants_count / data$Variants_count_sum
# Stacked
ggplot(data, aes(fill = Max_AF, y = Variants_count, x = Pathogenicity)) +
  geom_bar(position = "dodge", stat = "identity") + 
  scale_fill_nejm() +
  facet_grid(~ group) +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # Adjust font sizes
    axis.text.x = element_text(size = 14),     # Font size for x-axis labels
    axis.title.x = element_text(size = 16),    # Font size for x-axis title
    axis.text.y = element_text(size = 14),     # Font size for y-axis labels
    axis.title.y = element_text(size = 16),    # Font size for y-axis title
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  geom_text(aes(label = paste(Variants_count, " (", sprintf("%0.2f%%", Proportion * 100), ")", sep = "")), 
            position = position_dodge(width = 0.9), vjust = -0.5,size =4)

##Common plot  ===========
###Common clinvar variants plot

common_tab <- sorted_tab_old %>% filter(SZAID %in% intersect(unique(sorted_tab_old$SZAID),unique(swe_sorted_tab_old$SZAID)))
common_tab<- merge(common_tab,swe_sorted_tab_old[,c("SZAID","Freq")],by='SZAID',all.x= TRUE)
common_tab <- common_tab %>% filter(SZAID %in% intersect(unique(sorted_tab_old$SZAID),unique(swe_sorted_tab_old$SZAID)))
common_tab <- common_tab %>% select(SZAID, Freq.x,Freq.y,ClinVar_CLNSIG,Existing_variation,HGVSc,HGVSp,Database,MAX_AF)
common_tab$Freq.x <- common_tab$Freq.x/275
common_tab$Freq.y <- common_tab$Freq.y/101
# Assuming you have a dataframe called common_tab with columns Freq.x, Freq.y, and Subclass
# First, sort the dataframe by Freq.x in descending order and select the top 30 rows
top_30_features <- head(arrange(common_tab, desc(Freq.x)), 30)

ggplot(top_30_features, mapping = aes(x = SZAID, y = ClinVar_CLNSIG, fill = Freq.x)) +
  geom_tile() +
  geom_text(aes(label = Freq.x), vjust = 1) +  # Add labels for Freq.x values
  scale_fill_gradient(low = "white", high = "blue") +  # Adjust the color scale
  labs(x = "", y = "ClinVar_CLNSIG", fill = "Freq.x") +
  theme_minimal() +
  theme(axis.text.x = element_blank())  #

aggregated_data <- common_tab %>%
  mutate(MAX_AF_Category = case_when(
    is.na(MAX_AF) ~ "No_public AF",
    MAX_AF < 0.05 ~ "Public AF < 0.05",
    MAX_AF >= 0.05 ~ "Public AF >= 0.05"
  )) %>%
  group_by(MAX_AF_Category,ClinVar_CLNSIG) %>%
  summarize(Variants_count = n())
aggregated_data$Proportion <- aggregated_data$Variants_count / 180

aggregated_data <- aggregated_data%>%mutate(ClinVar_CLNSIG = 
         ifelse(ClinVar_CLNSIG == "Conflicting_interpretations_of_pathogenicity", "Conflicts", ClinVar_CLNSIG))

ggplot(aggregated_data, aes(fill = MAX_AF_Category, y = Variants_count, x = ClinVar_CLNSIG)) +
  geom_bar(position = "dodge", stat = "identity") + 
  scale_fill_manual(values = c("#EFC000FF","#868686FF"))  +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # Adjust font sizes
    axis.text.x = element_text(size = 14),     # Font size for x-axis labels
    axis.title.x = element_text(size = 16),    # Font size for x-axis title
    axis.text.y = element_text(size = 14),     # Font size for y-axis labels
    axis.title.y = element_text(size = 16),    # Font size for y-axis title
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  geom_text(aes(label = paste(Variants_count, " (", sprintf("%0.2f%%", Proportion * 100), ")", sep = "")), 
            position = position_dodge(width = 0.9), vjust = -0.5,size =4)

###both cohort most frequent ===========
Conflict <- sorted_tab %>% filter(ClinVar_CLNSIG=='Conflicting_interpretations_of_pathogenicity')%>% filter(Freq >13)
Conflict_swedish <- swe_sorted_tab_swedish %>% filter(ClinVar_CLNSIG=='Conflicting_interpretations_of_pathogenicity')%>% filter(Freq >5)
intersect(Conflict$SZAID,Conflict_swedish$SZAID)


AD_genes <- 
  unique(
    Summary_DB %>% filter(grepl("Autosomal dominant",inheritances))%>% 
      .$Genes
  )
AR_genes <- 
  unique(
    Summary_DB %>% filter(grepl("Autosomal recessive",inheritances))%>% 
      .$Genes
    )
##excluded the variants with conflict interpretations, any Clinvar variants that have public AF >= 0.05 together with mitochondrial variants-----------
sorted_tab_filter <- sorted_tab_old %>%filter(MAX_AF<0.05)%>%filter(ClinVar_CLNSIG!='Conflicting_interpretations_of_pathogenicity')
#sorted_tab_filter <- sorted_tab %>%filter(MAX_AF<0.05)%>%filter(ClinVar_CLNSIG!='Conflicting_interpretations_of_pathogenicity')
###这个不太对, MAX_AF should also included with NULL
sorted_tab_filter <-
  sorted_tab_filter %>%
  mutate(condition = case_when(Zygosity.x == "Homozygous" ~ "Positive",
                                               Genes %in% AD_genes & Zygosity.x == "Heterozygous" ~ "Positive",
                                               Genes %in% AR_genes & Zygosity.x == "Heterozygous" ~ "Carrier",
                                               TRUE ~ "Unsure"))
table(sorted_tab_filter[grepl('Basic',sorted_tab_filter$final_target_group),]$condition)
## old is the dataframe didn;t split the zygosity on each individual level count, not old one is the seoarated one
sorted_tab_old_filter <- sorted_tab_old %>%filter(MAX_AF<0.05)%>%filter(ClinVar_CLNSIG!='Conflicting_interpretations_of_pathogenicity')
# Group by Genes and calculate the sum of Freq for each gene
gene_counts <- sorted_tab_old_filter %>%
  group_by(Genes) %>%
  summarise(Total_Freq = sum(Freq))

# Arrange the result in descending order of Total_Freq
gene_counts <- gene_counts %>%
  arrange(desc(Total_Freq))

disease_counts <- sorted_tab_old_filter %>%
  group_by(Disease) %>%
  summarise(Total_Freq = sum(Freq))

# Print the table of genes with their total counts
print(gene_counts)
# Create a matrix of gene counts (rows = genes, columns = total frequency)
gene_matrix <- data.frame(Genes = gene_counts$Genes, Total_Freq = gene_counts$Total_Freq)

####venn diagram ========Figure 4B
x <- list(
  set1 <- unique(sorted_tab_filter$SZAID),
  set2 <- unique(sorted_tab_swedish_filter$SZAID)
)
display_venn(
  x,
  category.names = c('Turkish','Swedish'),
  fill = pal_jco("default")(2)
)

# Display the Venn diagram
grid.draw(venn.plot)

## ACMG 78 genes Healthy risk 
df_clinvar_filter <- df_clinvar %>% filter(MAX_AF<0.05)%>%filter(ClinVar_CLNSIG!='Conflicting_interpretations_of_pathogenicity')
df_clinvar_filter<-
  df_clinvar_filter %>%
  mutate(condition = case_when(Zygosity == "Homozygous" ~ "Positive",
                               Genes %in% AD_genes & Zygosity == "Heterozygous" ~ "Positive",
                               Genes %in% AR_genes & Zygosity== "Heterozygous" ~ "Carrier",
                               TRUE ~ "Unsure"))
#df_puta_filter <- df_puta %>% filter(MAX_AF<0.05)%>%filter(ClinVar_CLNSIG!='Conflicting_interpretations_of_pathogenicity')


ACMG <- df_clinvar_filter %>% filter(final_target_group=='Basic (for healthy subjects),Usually used for:Health predipositions/Disease risk')
write.xlsx(ACMG[,c('Genes','Disease','Existing_variation','HGVSc','Inheritances')],'/Users/xiyas/V2_Genome_reporting/table2.xlsx', sheetName = "Sheet1")


##===== putative variants
table(df_puta$MAX_AF < 0.05)

table(is.na(df_puta$MAX_AF))

####2: count the clinvar statictics among samples
puta_sorted_tab <- table(df_puta$Variant_info) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))

colnames(puta_sorted_tab) <- c('Variant_info','Freq')
test_puta <- df_puta[!duplicated(df_puta$Variant_info),]
test_puta$ClinVar_CLNSIG <- sapply(strsplit(test_puta$ClinVar_CLNSIG,"&"), `[`, 1)
puta_sorted_tab <- merge(x=puta_sorted_tab,y=test_puta,by='Variant_info',all.x = TRUE)
puta_sorted_tab <- puta_sorted_tab %>% arrange(desc(Freq))
#sorted_tab_puta$ClinVar_CLNSIG <- sapply(strsplit(sorted_tab_puta$ClinVar_CLNSIG,"&"), `[`, 1)
#sorted_tab$AC <- sapply(strsplit(sorted_tab$INFO,"&"), `[`, 1)

####remove mitochondrial genes ------------
puta_sorted_tab <- puta_sorted_tab %>% filter(X.CHROM != 'chrM')

####Category figure 5A
library(dplyr)
puta_sorted_tab <- puta_sorted_tab %>%
  mutate(group =
    case_when(IMPACT == 'HIGH' ~ "Impact_High(LoFs)",
              ada_score > 0.6 ~ "ada_score>0.6",
              rf_score > 0.6 ~ "rf_score >0.6",
              SpliceAI_pred_DS_AG > 0.5 ~ "SpliceAI_pred_DS_AG>0.5",
              SpliceAI_pred_DS_AL > 0.5 ~ "SpliceAI_pred_DS_AL>0.5",
              SpliceAI_pred_DS_DG > 0.5 ~ "SpliceAI_pred_DS_DG>0.5",
              SpliceAI_pred_DS_DL > 0.5 ~ "SpliceAI_pred_DS_DL>0.5",
              REVEL > 0.75 ~ "REVEL>0.75",
              BayesDel_addAF_score > 0.0692655 ~ "BayesDel_addAF_score >0.0692655 ",
              BayesDel_noAF_score > -0.0570105 ~ "BayesDel_noAF_score >-0.0570105"
    ),
    "Other"
  )

puta_sorted_tab <- puta_sorted_tab %>%
  mutate(MAX_AF_Category = case_when(
    is.na(MAX_AF) ~ "No_public AF",
    MAX_AF < 0.05 ~ "Public AF < 0.05",
    MAX_AF >= 0.05 ~ "Public AF >= 0.05"
  ))
puta_sorted_tab['Existing_variation'][puta_sorted_tab['Existing_variation'] ==''] <-NA


table(is.na(puta_sorted_tab$Existing_variation),puta_sorted_tab$MAX_AF_Category)
puta_sorted_tab <- puta_sorted_tab %>%
  mutate(Novel_or_Existing = case_when(
    is.na(Existing_variation) ~ "Novel variants",
    !is.na(Existing_variation) ~ "Existing variants"
  ))

puta_sorted_tab$group <- factor(puta_sorted_tab$group)
novel_conse <- puta_sorted_tab %>% group_by(Novel_or_Existing, Consequence) %>%
  summarize(count = n())
# Reorder the levels of the "group" factor variable
puta_sorted_tab$group <- relevel(puta_sorted_tab$group, ref = "Impact_High(LoFs)")

#Get the contingency table
contingency_table <- table(puta_sorted_tab$Novel_or_Existing,puta_sorted_tab$group, puta_sorted_tab$MAX_AF_Category)

# Print the contingency table
print(contingency_table)
contingency_df <- as.data.frame(contingency_table)

# Print the data frame
print(contingency_df)
# Create a scatter plot with three variables
ggplot(contingency_df, aes(x = Var3, y = Var2, size=Freq,color = Var1)) +
  geom_point() +
  geom_text(aes(label = Freq,size=200), vjust = -1.5) +
  scale_color_manual(values = c("royalblue3", "indianred")) +  # Set custom colors for Z levels
  labs(x = "X-Axis", y = "Y-Axis", color = "Z-Axis") +  # Set axis labels
  facet_grid(. ~Var1) +
  xlab("public AF")+ylab("in silico group")+
  theme_minimal()
  

puta_sorted_tab$CLIN_SIG <- sapply(strsplit(puta_sorted_tab$CLIN_SIG,"/"), `[`, 1)
puta_sorted_tab$CLIN_SIG <- sapply(strsplit(puta_sorted_tab$CLIN_SIG,"&"), `[`, 1)
puta_sorted_tab <- puta_sorted_tab %>% mutate(condition = case_when(Zygosity == "Homozygous" ~ "Positive",
                             Genes %in% AD_genes & Zygosity == "Heterozygous" ~ "Positive",
                             Genes %in% AR_genes & Zygosity == "Heterozygous" ~ "Carrier",
                             TRUE ~ "Negative"))
impact_High <- puta_sorted_tab%>%filter(group == "Impact_High(LoFs)")
##singletons? 
table(impact_High$Freq ==1 &impact_High$Zygosity =='Heterozygous')
impact_High['Existing_variation'][impact_High['Existing_variation'] ==''] <-NA

table(is.na(impact_High$Existing_variation))
#singleton
table(is.na(impact_High$Existing_variation) & impact_High$Freq ==1 &impact_High$Zygosity =='Heterozygous' )
table(impact_High$Zygosity =='Homozygous')
table(is.na(impact_High$Existing_variation) & impact_High$Zygosity =='Homozygous')

##ACMG =======
table(impact_High$final_target_group == 'Basic (for healthy subjects),Usually used for:Health predipositions/Disease risk')
#table_2 <- impact_High[impact_High$final_target_group == 'Basic (for healthy subjects),Usually used for:Health predipositions/Disease risk',]
ACMG_LoFs <- impact_High%>% 
  filter(final_target_group == 'Basic (for healthy subjects),Usually used for:Health predipositions/Disease risk')
table(is.na(ACMG_LoFs$Existing_variation))
table(ACMG_LoFs$condition)
ada_score
dim(impact_High)


# check overlaps
group_variants <- list()
# Loop through each group category
for (group_category in unique(puta_sorted_tab$group)) {
  group_variants[[group_category]] <- puta_sorted_tab %>%
    filter(group == group_category) %>%
    select(Variant_info) %>%
    pull()
}
overlapping_variants <- Reduce(intersect, group_variants)

# If you want to see the overlapping variants, you can print them:
print(overlapping_variants)


############grouped-pie-chart ============
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
  AD_genes <- 
    unique(
      Summary_DB %>% filter(grepl(KEY,Target.group)) %>% 
        filter(grepl("Autosomal dominant",inheritances))%>% 
        .$Genes
    )
  AR_genes <- 
    unique(
      Summary_DB %>% filter(grepl(KEY,Target.group)) %>% 
        filter(grepl("Autosomal recessive",inheritances))%>% 
        .$Genes
    )
  aa <- 
    df_clinvar_filter%>% filter(grepl(KEY,final_target_group)) %>%
    mutate(condition = case_when(Zygosity == "Homozygous" ~ "Positive",
                                 Genes %in% AD_genes & Zygosity == "Heterozygous" ~ "Positive",
                                 Genes %in% AR_genes & Zygosity == "Heterozygous" ~ "Carrier",
                                 TRUE ~ "Negative"))
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
  #pie(c(count,100-count),labels = c("Carrier(n=39)","Non-Carrier(n=61)"),col = c("#F6AE2D", "#33658A"))
}
#########based on putative variants then too much --------------
bb <- df_puta%>% filter(Genes %in% target_genes) %>%
  mutate(condition = case_when(Genes %in% AD_genes ~ "Positive",
                               Genes %in% AR_genes | Zygosity == "Homozygous" ~ "Positive",
                               TRUE ~ "Negative"))

unique(bb %>% filter(condition == "Positive") %>% .$patientID)

####I just want to count the positive results patients vs total patients, based on clinvar variants
##This I try with pie chart of different panels ------------------
rm(test)
library(ggplot2)

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

########figure of grouped bar plot --------
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

###start count the Positive Negative and Carrier -----------
###should not use genes, should use groups
AD_genes <- 
  unique(
    Summary_DB %>% filter(grepl("Autosomal dominant",inheritances))%>% 
      .$Genes
  )
AR_genes <- 
  unique(
    Summary_DB %>% filter(grepl("Autosomal recessive",inheritances))%>% 
      .$Genes
  )
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
  scale_y_continuous(expand = c(0,0), limits = c(0, 6)) +
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


####Case report =========
P001_106 <-df_clinvar_filter[df_clinvar_filter$patientID == "P001_106.hard-filtered.vcf.gz_vep_annotated.vcf_outFile_sp_Inheritance_4_nodup.txt",]
P001_85 <- df_clinvar_filter[df_clinvar_filter$patientID == "P001_85.hard-filtered.vcf.gz_vep_annotated.vcf_outFile_sp_Inheritance_4_nodup.txt",]
P001_96 <- df_clinvar_filter[df_clinvar_filter$patientID == "P001_96.hard-filtered.vcf.gz_vep_annotated.vcf_outFile_sp_Inheritance_4_nodup.txt",]

P001_106_more <-df_clinvar[df_clinvar$patientID == "P001_106.hard-filtered.vcf.gz_vep_annotated.vcf_outFile_sp_Inheritance_4_nodup.txt",]
P001_85_more <- df_clinvar[df_clinvar$patientID == "P001_85.hard-filtered.vcf.gz_vep_annotated.vcf_outFile_sp_Inheritance_4_nodup.txt",]
P001_96_more <- df_clinvar[df_clinvar$patientID == "P001_96.hard-filtered.vcf.gz_vep_annotated.vcf_outFile_sp_Inheritance_4_nodup.txt",]

