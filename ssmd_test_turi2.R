#trying the SSMD method on turi dataframe
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics")

#load packages that I'll probably need
library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(lspline)

#for DMSO only
turi2_df_bac_DMSO_spread$diff <- turi2_df_bac_DMSO_spread$"8" -  turi2_df_bac_DMSO_spread$"1"
turi2_df_bac_DMSO_avg <- mean(turi2_df_bac_DMSO_spread$diff)
turi2_df_library_DMSO_spread$diff <- turi2_df_library_DMSO_spread$"8" -  turi2_df_library_DMSO_spread$"1"
turi2_df_library_DMSO_spread$diff_neg <- turi2_df_library_DMSO_spread$diff - turi2_df_bac_DMSO_avg
turi2_df_library_DMSO_spread_bywells <- group_by(turi2_df_library_DMSO_spread, well) %>% summarize(m = mean(diff_neg))
turi2_df_library_DMSO_spread_bywells_test <- left_join(turi2_df_library_DMSO_spread, turi2_df_library_DMSO_spread_bywells, by = "well")    
turi2_df_library_DMSO_spread_bywells_test$sub <- turi2_df_library_DMSO_spread_bywells_test$diff_neg - turi2_df_library_DMSO_spread_bywells_test$m
turi2_df_library_DMSO_spread_bywells_test$sub2 <- (turi2_df_library_DMSO_spread_bywells_test$sub)^2
turi2_df_library_DMSO_spread_scalc <- group_by(turi2_df_library_DMSO_spread_bywells_test, well) %>% summarize(sum = sum(sub2))
turi2_df_library_DMSO_spread_scalc$s_val <- sqrt((turi2_df_library_DMSO_spread_scalc$sum)/2)
turi2_df_library_DMSO_spread_scalc_dcalc <- left_join(turi2_df_library_DMSO_spread_scalc, turi2_df_library_DMSO_spread_bywells_test, by = "well")    
turi2_df_library_DMSO_spread_scalc_dcalc$mm <- turi2_df_library_DMSO_spread_scalc_dcalc$m/turi2_df_library_DMSO_spread_scalc_dcalc$s_val
#need to determine cutoff