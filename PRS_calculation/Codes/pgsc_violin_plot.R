library(ggplot2)
library(RColorBrewer)
# Load data from txt files

data2 <- read.table("/Users/xiyas/pgsc_calc/Turkish_4_score/results/score/aggregated_scores.txt.gz",header = TRUE)
data3 <- read.table("/Users/xiyas/pgsc_calc/Swedish_4_score/results/score/aggregated_scores.txt.gz",header = TRUE)



palette <- brewer.pal(n = 6, name = "Pastel1")[1:3]
# Combine data into one data frame
combined_data <- data.frame(value = c(data2$PGS002758_hmPOS_GRCh38_AVG, data3$PGS002758_hmPOS_GRCh38_AVG), 
                            group = factor(rep(c("Turkish", "Swedish"), 
                                               times = c(nrow(data2), nrow(data3)))))

# Create violin plot
ggplot(combined_data, aes(x = group, y = value,fill=group)) + 
  geom_violin() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  scale_fill_manual(values = palette) +
  labs(x = "Data Group", y = "Value", title = "PGS000014_AVG Distribution of two Datasets") +
  theme_bw()

