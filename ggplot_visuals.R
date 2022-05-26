#sanity checking

#run necessary packages
library(tidyverse)
library(dplyr)
library(ggplot2)

#can look at the combination of all the graphs with unique "id" value
all_curves_solvent <- function(df){
  df_solv_gc <- spread(df, key = time, value = OD)
  col_num <- ncol(df_solv_gc)
  df_solv_gc <- select(df_solv_gc, well, rep, c(6:col_num))
  df_solv_gc$id <- paste(df_solv_gc$well, df_solv_gc$rep)
  df_solv_gc <- select(df_solv_gc, c(3:12))
  df_solv_gc <- melt(df_solv_gc, id.vars = 'id')
  
  
  return(ggplot(df_solv_gc, aes(x = variable, y = value)) + 
    geom_line(aes(color = id, group = id)) + 
    theme(legend.position="none"))
}
##############################3

#look at specific things of interest

###rep1
rep1_solvent <- function(df){
  df_solv_gc <- spread(df, key = time, value = OD)
  col_num <- ncol(df_solv_gc)
  df_solv_gc <- select(df_solv_gc, well, rep, c(6:col_num))
  df_solv_gc$id <- paste(df_solv_gc$well, df_solv_gc$rep)
  df_solv_gc <- select(df_solv_gc, c(3:12))
  df_solv_gc <- melt(df_solv_gc, id.vars = 'id')
  df_solv_gc <-  df_solv_gc[id %like% "rep1"]
  
  
  return(ggplot(df_solv_gc, aes(x = variable, y = value)) + 
           geom_line(aes(color = id, group = id)) + 
           theme(legend.position="none"))
}
##rep2
rep2_solvent <- function(df){
  df_solv_gc <- spread(df, key = time, value = OD)
  col_num <- ncol(df_solv_gc)
  df_solv_gc <- select(df_solv_gc, well, rep, c(6:col_num))
  df_solv_gc$id <- paste(df_solv_gc$well, df_solv_gc$rep)
  df_solv_gc <- select(df_solv_gc, c(3:12))
  df_solv_gc <- melt(df_solv_gc, id.vars = 'id')
  df_solv_gc <-  df_solv_gc[id %like% "rep2"]
  
  
  return(ggplot(df_solv_gc, aes(x = variable, y = value)) + 
           geom_line(aes(color = id, group = id)) + 
           theme(legend.position="none"))
}
##rep3 
rep3_solvent <- function(df){
  df_solv_gc <- spread(df, key = time, value = OD)
  col_num <- ncol(df_solv_gc)
  df_solv_gc <- select(df_solv_gc, well, rep, c(6:col_num))
  df_solv_gc$id <- paste(df_solv_gc$well, df_solv_gc$rep)
  df_solv_gc <- select(df_solv_gc, c(3:12))
  df_solv_gc <- melt(df_solv_gc, id.vars = 'id')
  df_solv_gc <-  df_solv_gc[id %like% "rep3"]
  
  
  return(ggplot(df_solv_gc, aes(x = variable, y = value)) + 
           geom_line(aes(color = id, group = id)) + 
           theme(legend.position="none"))
}
##singles note: need to pass the well name you're interested in in quotes, and it'll work
singles_solvent <- function(df, well){
  df_solv_gc <- spread(df, key = time, value = OD)
  col_num <- ncol(df_solv_gc)
  df_solv_gc <- select(df_solv_gc, well, rep, c(6:col_num))
  df_solv_gc$id <- paste(df_solv_gc$well, df_solv_gc$rep)
  df_solv_gc <- select(df_solv_gc, c(3:12))
  df_solv_gc <- melt(df_solv_gc, id.vars = 'id')
  df_solv_gc <-  df_solv_gc[id %like% {{well}}]
  
  
  return(ggplot(df_solv_gc, aes(x = variable, y = value)) + 
           geom_line(aes(color = id, group = id)) + 
           theme(legend.position="none"))
}

