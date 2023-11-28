read_files <- function(path) {
  setwd(path)
  files <- list.files(path = path, pattern = "_pharmaco.txt")
  list_temp <- lapply(files, function(file) {
    temp <- read.delim(file)
    temp <- temp[!duplicated(temp$SZAvarID), ]
    temp$patientID <- rep(file)
    return(temp)
  })
  return(do.call(rbind, list_temp))
}

df_pharmaco_TR <- read_files("/Users/xiyas/V2_Genome_reporting/PGx_results_TR") 
df_pharmaco_SW <- read_files("/Users/xiyas/V2_Genome_reporting/PGx_results_SW") 
#df_pharmaco_TR %>% group_by(patientID) %>% filter(Level.Modifiers == 'Tier 1 VIP') %>% summarise(Count = n())
count <- df_pharmaco_TR %>% group_by(patientID) %>% filter(Level.Modifiers == 'Tier 1 VIP') %>%
  filter(Level.of.Evidence %in% c('1A', '1B', '2A','3')) %>% summarise(Count = n())
median <- median(count$Count)
median
count <- df_pharmaco_TR %>% group_by(patientID) %>% filter(Level.Modifiers == 'Tier 1 VIP') %>%
  filter(Level.of.Evidence %in% c('1A', '1B', '2A')) %>% summarise(Count = n())

count <- df_pharmaco_SW %>% group_by(patientID) %>% filter(Level.Modifiers == 'Tier 1 VIP') %>%
  filter(Level.of.Evidence %in% c('1A', '1B', '2A','3')) %>% summarise(Count = n())
median <- median(count$Count)
median
count <- df_pharmaco_SW %>% group_by(patientID) %>% filter(Level.Modifiers == 'Tier 1 VIP') %>%
  filter(Level.of.Evidence %in% c('1A', '1B', '2A')) %>% summarise(Count = n())
median <- median(count$Count)
median
