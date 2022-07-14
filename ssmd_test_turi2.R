#trying the SSMD method on turi dataframe
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics")

#load packages that I'll probably need
library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(lspline)
library("fuzzyjoin")


#for DMSO only
#bascially, calculate the values needed to find the SSMD (see XD Zhang paper for explanation)
#I calculate the MM SSMD here
bac$diff <- bac$"8" -  bac$"1"
turi2_df_bac_DMSO_avg <- mean(bac$diff)
lib$diff <- lib$"8" -  lib$"1"
lib$diff_neg <- lib$diff - turi2_df_bac_DMSO_avg
lib_bywells <- group_by(lib, well) %>% summarize(m = mean(diff_neg))
lib_bywells_test <- left_join(lib, lib_bywells, by = "well")    
lib_bywells_test$sub <- lib_bywells_test$diff_neg - lib_bywells_test$m
lib_bywells_test$sub2 <- (lib_bywells_test$sub)^2
lib_scalc <- group_by(lib_bywells_test, well) %>% summarize(sum = sum(sub2))
lib_scalc$s_val <- sqrt((lib_scalc$sum)/2)
lib_scalc_dcalc <- left_join(lib_scalc, lib_bywells_test, by = "well")    
lib_scalc_dcalc$mm <- lib_scalc_dcalc$m/lib_scalc_dcalc$s_val
#need to determine cutoff

#write function for other datasets?
#need to input the spread and filtered version of bac/library dfs -> the name of the dataframe should NOT contain "spread"
ssmd_test <- function(bac, lib){
  bac$diff <- bac$"8" -  bac$"1"
  turi2_df_bac_DMSO_avg <- mean(bac$diff)
  lib$diff <- lib$"8" -  lib$"1"
  lib$diff_neg <- lib$diff - turi2_df_bac_DMSO_avg
  lib_bywells <- group_by(lib, well) %>% summarize(m = mean(diff_neg))
  lib_bywells_test <- left_join(lib, lib_bywells, by = "well")    
  lib_bywells_test$sub <- lib_bywells_test$diff_neg - lib_bywells_test$m
  lib_bywells_test$sub2 <- (lib_bywells_test$sub)^2
  lib_scalc <- group_by(lib_bywells_test, well) %>% summarize(sum = sum(sub2))
  lib_scalc$s_val <- sqrt((lib_scalc$sum)/2)
  lib_scalc_dcalc <- left_join(lib_scalc, lib_bywells_test, by = "well")    
  lib_scalc_dcalc$mm <- lib_scalc_dcalc$m/lib_scalc_dcalc$s_val
  lib_scalc_dcalc <- lib_scalc_dcalc[!duplicated(lib_scalc_dcalc$well),]
  return(lib_scalc_dcalc)
}

#turns out that using the filtered version for turi does make a difference on the results, since there are some wells without growth
#test_turi2_MeOH_2 <- ssmd_test(turi2_df_bac_MeOH_spread, turi2_df_library_MeOH_spread)
test_turi2_MeOH <- ssmd_test(turi2_df_bac_MeOH, turi2_df_library_MeOH)
test_turi2_H2O <- ssmd_test(turi2_df_bac_H2O, turi2_df_library_H2O)
test_turi2_DMSO <- ssmd_test(turi2_df_bac_DMSO, turi2_df_library_DMSO)

#join with platemap?
test_turi2_H2O_map <- left_join(test_turi2_H2O, platemap_H2O_merge_form, by = "well")
#we need the chemical names, so use the manually edited file that essentially maps the messed up names with the regex names
# THERE HAS TO BE A BETTER WAY OF DOING THIS, LIKE ASSIGNING A CODE OR SMTH, BUT IM DUMB and have no foresight haha
turi_jw_key <- read.csv("turi_out_fuzzy_jw_edited_DO_NOT_EDIT.csv")
colnames(turi_jw_key) <- c("X", "solv.x", "well", "category", "solv", "dist")
test_turi2_H2O_map <- left_join(test_turi2_H2O_map, turi_jw_key, by = "solv") 
#assign color based on mm value
colfunc <- colorRampPalette(c('purple', 'green'))

#essentially using the order of the the "mm" value, arrange and assign a rank column with the row number
#testing for function########
test_turi2_H2O_map_color <- test_turi2_H2O_map %>%
  drop_na(mm) %>% 
  arrange(mm) %>%
  mutate(rank = row_number())

#join the initial dataframe and the dataframe with the color gradient by the mm value, and "rank" the colors generated 
#i'm not 100% sure how this works either
test_turi2_H2O_map <- left_join(test_turi2_H2O_map, test_turi2_H2O_map_color, by = "mm")
the_colors <- tibble(rank = 1:max(test_turi2_H2O_map$rank, na.rm = T), color_code = colfunc(max(test_turi2_H2O_map$rank, na.rm = T)))

test_turi2_H2O_map_color_test <- test_turi2_H2O_map %>%
  left_join(the_colors, by = 'rank') %>%
  mutate(color_code = if_else(is.na(mm), '#FFFFFFFF', color_code))

#select out the columns of interest from the mess of columns that is a result of joining over and over again
test_turi2_H2O_map_color_test <- test_turi2_H2O_map_color_test[!duplicated(test_turi2_H2O_map_color_test$solv.x.y),]
test_turi2_H2O_map_color_test_clean <- select(test_turi2_H2O_map_color_test, solv.x.y, color_code)

#make into function#########
#need to make a separately formatted sheet for H2O, since there are some special characters that don't get imported correctly
#there are mismatches otherwise, but just put this platemap into the function instead of the one generated earlier
platemap_H2O_merge_form <- read.csv("platemap_H2O_merge.csv")

color_function <- function(df, map){
  df_map <- left_join(df, map, by = "well") 
  turi_jw_key <- read.csv("turi_out_fuzzy_jw_edited_DO_NOT_EDIT.csv")
  colnames(turi_jw_key) <- c("X", "solv", "well", "category", "solv", "dist")
  df_map <- left_join(df_map, turi_jw_key, by = "solv") 
  #assign color based on mm value
  colfunc <- colorRampPalette(c('purple', 'green'))
  
  df_map_color <- df_map %>%
    drop_na(mm) %>% 
    arrange(mm) %>%
    mutate(rank = row_number())
  
  df_map <- left_join(df_map, df_map_color, by = 'mm')
  the_colors <- tibble(rank = 1:max(df_map$rank, na.rm = T), color_code = colfunc(max(df_map$rank, na.rm = T)))
  
  df_map_color_test <- df_map %>%
    left_join(the_colors, by = 'rank') %>%
    mutate(color_code = if_else(is.na(mm), '#FFFFFFFF', color_code))
  df_map_color_test <- df_map_color_test[!duplicated(df_map_color_test$solv.x.y),]
  df_map_color_test_clean <- select(df_map_color_test, solv.x.y, color_code, well.y.y)
  return(df_map_color_test_clean)
}

#output each color assigned to each chemical
turi2_DMSO_color_out <- color_function(test_turi2_DMSO, platemap_DMSO_merge)
turi2_MeOH_color_out <- color_function(test_turi2_MeOH, platemap_MeOH_merge)
turi2_H2O_color_out <- color_function(test_turi2_H2O, platemap_H2O_merge_form)
#combine into one df
turi2_color_out <- rbind(turi2_DMSO_color_out, turi2_MeOH_color_out, turi2_H2O_color_out)

#match column names so that 
colnames(turi2_color_out) <- c("Chemical", "color_code", "well")

#import the itol file so the names used by iTOL can be assigned to the names from regex
ssmd_turi_itol <- read.csv("iTOL_files/dataset_symbols_template_colors_test_copy.txt")

#have to use fuzzy join since some characters are slightly added or removed
turi2_color_out_labels <- stringdist_join(turi2_color_out, ssmd_turi_itol, method = "jw", by = "Chemical", max_dist = 99, distance_col = "dist") %>% group_by(Chemical.x) %>% slice_min(order_by = dist, n=1)
turi2_color_out_labels <- data.frame(turi2_color_out_labels)
turi2_color_out_labels <- select(turi2_color_out_labels, Chemical.y, Symbol, Size, color_code, Fill, Position)

#output file is very close to iTOL useable, just need to paste into the template
write.csv(turi2_color_out_labels, "ssmd_turi2_color_out.csv", quote = F, row.names = F)

#manually setting the color values?
test_turi2_DMSO_color_2 <- test_turi2_DMSO
test_turi2_MeOH_color_2 <- test_turi2_MeOH
test_turi2_H2O_color_2 <- test_turi2_H2O

#using nested ifelse, find values above and below the cutoffs\
#there is probably a better way to do this using "when"
####note to self I think there is an extra step in here
test_turi2_DMSO_color_2$color_2 <- ifelse(test_turi2_DMSO_color_2$mm <= -1.645, "#8b12b3", 
                                          ifelse(test_turi2_DMSO_color_2$mm >=1.645, "#12b31d",
                                                 ifelse(test_turi2_DMSO_color_2$mm > -1.645 & test_turi2_DMSO_color_2$mm < 1.645, "#8f8f8f", "NA")))
test_turi2_H2O_color_2$color_2 <- ifelse(test_turi2_H2O_color_2$mm <= -1.645, "#8b12b3", 
                                         ifelse(test_turi2_H2O_color_2$mm >=1.645, "#12b31d",
                                                ifelse(test_turi2_H2O_color_2$mm > -1.645 & test_turi2_H2O_color_2$mm < 1.645, "#8f8f8f", "NA")))
test_turi2_MeOH_color_2$color_2 <- ifelse(test_turi2_MeOH_color_2$mm <= -1.645, "#8b12b3", 
                                          ifelse(test_turi2_MeOH_color_2$mm >=1.645, "#12b31d",
                                                 ifelse(test_turi2_MeOH_color_2$mm > -1.645 & test_turi2_MeOH_color_2$mm < 1.645, "#8f8f8f", "NA")))


#correlate to chemical names from regex
test_turi2_DMSO_color_2 <- left_join(test_turi2_DMSO_color_2, platemap_DMSO_merge, by = "well")
test_turi2_DMSO_color_2 <- select(test_turi2_DMSO_color_2, well, mm, color_2, solv.x)
colnames(test_turi2_DMSO_color_2) <- c("well", "mm", "color_2", "solv")

test_turi2_MeOH_color_2 <- left_join(test_turi2_MeOH_color_2, platemap_MeOH_merge, by = "well")
test_turi2_MeOH_color_2 <- select(test_turi2_MeOH_color_2, well, mm, color_2, solv)

test_turi2_H2O_color_2 <- left_join(test_turi2_H2O_color_2, platemap_H2O_merge_form, by = "well")
test_turi2_H2O_color_2 <- select(test_turi2_H2O_color_2, well, mm, color_2, solv)

#combine
test_turi2_color_2 <- rbind(test_turi2_DMSO_color_2, test_turi2_H2O_color_2, test_turi2_MeOH_color_2)
colnames(test_turi2_color_2) <- c("well", "mm", "color_2", "solv")


colnames(turi_jw_key) <- c("X", "solv.2", "well", "category", "solv", "dist")
test_turi2_color_2_jw <- stringdist_join(test_turi2_color_2, turi_jw_key, method = "jw", by = "solv", max_dist = 99, distance_col = "dist_2") %>% group_by(solv.x) %>% slice_min(order_by = dist_2, n=1)

test_turi2_color_2_jw <- data.frame(test_turi2_color_2_jw)
test_turi2_color_2 <- select(test_turi2_color_2_jw, well.x, color_2, solv.2)
colnames(test_turi2_color_2) <- c("well", "color_2", "Chemical")

#coordinate to the names used in iTOL
turi2_color_out_labels_2 <- stringdist_join(test_turi2_color_2, ssmd_turi_itol, method = "jw", by = "Chemical", max_dist = 99, distance_col = "dist") %>% group_by(Chemical.x) %>% slice_min(order_by = dist, n=1)
turi2_color_out_labels_2_out <- select(data.frame(turi2_color_out_labels_2), Chemical.y, Symbol, Size, color_2, Fill, Position)
write.csv(turi2_color_out_labels_2_out, "ssmd_turi2_color_out_2.csv", quote = F, row.names = F)
