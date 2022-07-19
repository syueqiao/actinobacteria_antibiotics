library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(lspline) 

setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/media_test/TochevaData/cyto12")
cyto_num <- list.files(path="./", full.names = T , recursive =F)
method <- list.files(path=paste0(cyto_num), full.names = T , recursive =F)
files <- list.files(path=paste0(method), pattern="*.CSV", full.names=TRUE, recursive=FALSE)
datalist = list()

#make dataframe that contains all the files with metadata appended

df <- data.frame()
for(file in files){
  new_df <- fread(file, select = c(1:3))
  new_df$id <- rep(file, nrow(new_df))
  df <- rbind(df, new_df, use.names=FALSE)
}

df %>% separate(id, c(NA ,"bug", "solvent", "rep", "filename"), sep = "([/])") -> df_metadata
df_metadata %>% separate(filename, c(NA, NA, "t"), sep = "([_])") -> df_metadata
df_metadata$t <- gsub('.{4}$', '', df_metadata$t)
df_metadata$t <- sub('.', '', df_metadata$t)
colnames(df_metadata) <- c("well", "content", "OD", "bug", "solvent", "rep", "time") 

