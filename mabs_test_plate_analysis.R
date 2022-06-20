#load packages
library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(lspline)

setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/mabs_test")

#import platemap
m_platemaps <- platemap <- read.csv("../m_platemaps.csv", header = TRUE)

#import the data
files <- list.files(path="./", full.names = T , recursive =F)

#makes all the data into one dataframe, with time listed in the columns
df <- data.frame()
for(file in files){
  new_df <- fread(file, select = c(1:3))
  new_df$id <- rep(file, nrow(new_df))
  df <- rbind(df, new_df, use.names=FALSE)
}

df %>% separate(id, c(NA, NA, "t"), sep = "([_])") -> df_metadata
df_metadata$t <- gsub('.{4}$', '', df_metadata$t)
df_metadata$t <- sub('.', '', df_metadata$t)
colnames(df_metadata) <- c("well", "content", "OD", "time")

#set lists for each of the solvents/contains/growth
#FOR THE MEDIA#####
mabs_media_list <- unique(m_platemaps$Growth)

media_out <- list()
for (media in mabs_media_list) {
 filtered <- filter(m_platemaps, Growth == media)
 media_out[[media]] <- filtered
  
}

#retrieve each list
media_7TOGLY <- media_out$`7TOGLY`$Well
media_LB <- media_out$`LB`$Well
media_7TGLY <- media_out$`7TGLY`$Well
media_7TO <- media_out$`7TO`$Well
media_LTOGLU <- media_out$`LTOGLU`$Well
media_LTO <- media_out$`LTO`$Well
media_LT <- media_out$`LT`$Well

#FOR THE SOLVENTS#####
mabs_solv_list <- unique(m_platemaps$Solvent)

solv_out <- list()
for (solv in mabs_solv_list) {
  filtered <- filter(m_platemaps, Solvent == solv)
  solv_out[[solv]] <- filtered
  
}

#retrieve each list
solv_none <- solv_out$`none`$Well
solv_DMSO <- solv_out$`DMSO`$Well
solv_H2O <- solv_out$`H2O`$Well
solv_MeOH <- solv_out$`MeOH`$Well
solv_DMF <- solv_out$`DMF`$Well
solv_EtOH <- solv_out$`EtOH`$Well

#FOR WELL CONTENTS#######
mabs_content_list <- unique(m_platemaps$Contents)

content_out <- list()
for (cont in mabs_content_list) {
  filtered <- filter(m_platemaps, Contents == cont)
  content_out[[cont]] <- filtered
  
}

content_bac <- content_out$`Bacteria`$Well
content_kana <- content_out$`Kana`$Well
content_emp <- content_out$`Empty`$Well

#make a vector for when you don't want to filter by anything in a certain category
all <- m_platemaps$Well

#visualize growth curves  using ggplot
#something like this, may have to be spread

mabs_curves <- function(df, c, m, so){
  
  df_in <- filter(df, well %in% c)
  df_in <- filter(df_in, well %in% m)
  df_in <- filter(df_in, well %in% so)
  
    
    
return(ggplot(df_in, aes(x = time, y = OD)) + 
         geom_line(aes(color = well, group =well)) +
         ylim(0, 1) +
         ggtitle(substitute(c), substitute(m)) +
         theme(legend.position="none"))
}
#plot together, each well should have 6 replicates for each solv/media combination
#first input is dataframe, next is the content of the wells, finall the solvent

mabs_curves(df_metadata, content_bac, media_7TOGLY, solv_DMSO)
