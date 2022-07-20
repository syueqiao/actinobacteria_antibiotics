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
#bac section
###
new_ssmd <- function(bac_in, lib_in){
bac_in$diff <- bac_in$"8" -  bac_in$"1"

turi2_df_bac_in_DMSO_avg_rep1 <- bac_in %>%
  filter(rep == "rep1") %>%
  .$diff %>%
  mean

turi2_df_bac_in_DMSO_avg_rep2 <- bac_in %>%
  filter(rep == "rep2") %>%
  .$diff %>%
  mean

turi2_df_bac_in_DMSO_avg_rep3 <- bac_in %>%
  filter(rep == "rep3") %>%
  .$diff %>%
  mean

lib_in$diff <- lib_in$"8" -  lib_in$"1"

turi2_df_lib_in_DMSO_avg_rep1 <- lib_in %>%
  filter(rep == "rep1")
turi2_df_lib_in_DMSO_avg_rep1$d_1 <- turi2_df_lib_in_DMSO_avg_rep1$diff -  turi2_df_bac_in_DMSO_avg_rep1

turi2_df_lib_in_DMSO_avg_rep2 <- lib_in %>%
  filter(rep == "rep2")
turi2_df_lib_in_DMSO_avg_rep2$d_1 <- turi2_df_lib_in_DMSO_avg_rep2$diff -  turi2_df_bac_in_DMSO_avg_rep2

turi2_df_lib_in_DMSO_avg_rep3 <- lib_in %>%
  filter(rep == "rep3")
turi2_df_lib_in_DMSO_avg_rep3$d_1 <- turi2_df_lib_in_DMSO_avg_rep3$diff -  turi2_df_bac_in_DMSO_avg_rep3

turi2_df_lib_in_DMSO_avg_repall <- rbind (turi2_df_lib_in_DMSO_avg_rep1, turi2_df_lib_in_DMSO_avg_rep2, turi2_df_lib_in_DMSO_avg_rep3)
turi2_df_lib_in_DMSO_avg_repall_avg <- group_by(turi2_df_lib_in_DMSO_avg_repall, well) %>%
  summarize(d_calc = mean(d_1)) %>%
  left_join(., turi2_df_lib_in_DMSO_avg_repall, by = "well")

turi2_df_lib_in_DMSO_avg_repall_avg$s_1 <- (turi2_df_lib_in_DMSO_avg_repall_avg$d_1 - turi2_df_lib_in_DMSO_avg_repall_avg$d_calc)^2

turi2_df_lib_in_DMSO_avg_repall_s_2 <- group_by(turi2_df_lib_in_DMSO_avg_repall_avg, well)%>%
  summarize(s_3 = sum(s_1))

turi2_df_lib_in_DMSO_avg_repall_s_2$s_calc <- sqrt(turi2_df_lib_in_DMSO_avg_repall_s_2$s_3/2)
turi2_df_lib_in_DMSO_avg_repall_avg <- left_join(turi2_df_lib_in_DMSO_avg_repall_s_2, turi2_df_lib_in_DMSO_avg_repall_avg, by = "well")
turi2_df_lib_in_DMSO_avg_repall_avg$mm <- turi2_df_lib_in_DMSO_avg_repall_avg$d_calc/turi2_df_lib_in_DMSO_avg_repall_avg$s_calc
turi2_df_lib_in_DMSO_avg_repall_avg <- turi2_df_lib_in_DMSO_avg_repall_avg[!duplicated(turi2_df_lib_in_DMSO_avg_repall_avg$well),]
return(turi2_df_lib_in_DMSO_avg_repall_avg)
}

turi2_DMSO_ssmd <- new_ssmd(turi2_df_bac_DMSO, turi2_df_library_DMSO)
turi2_MeOH_ssmd <- new_ssmd(turi2_df_bac_MeOH, turi2_df_library_MeOH)
turi2_H2O_ssmd <- new_ssmd(turi2_df_bac_H2O, turi2_df_library_H2O)

##################################################

#turi2_df_library_DMSO$diff_neg <- turi2_df_library_DMSO$diff - turi2_df_turi2_df_bac_DMSO_DMSO_avg
#turi2_df_library_DMSO_bywells <- group_by(turi2_df_library_DMSO, well) %>% summarize(m = mean(diff_neg))
#turi2_df_library_DMSO_bywells_test <- left_join(turi2_df_library_DMSO, turi2_df_library_DMSO_bywells, by = "well")    
#turi2_df_library_DMSO_bywells_test$sub <- turi2_df_library_DMSO_bywells_test$diff_neg - turi2_df_library_DMSO_bywells_test$m
#turi2_df_library_DMSO_bywells_test$sub2 <- (turi2_df_library_DMSO_bywells_test$sub)^2
#turi2_df_library_DMSO_scalc <- group_by(turi2_df_library_DMSO_bywells_test, well) %>% summarize(sum = sum(sub2))
#turi2_df_library_DMSO_scalc$s_val <- sqrt((turi2_df_library_DMSO_scalc$sum)/2)
#turi2_df_library_DMSO_scalc_dcalc <- left_join(turi2_df_library_DMSO_scalc, turi2_df_library_DMSO_bywells_test, by = "well")    
#turi2_df_library_DMSO_scalc_dcalc$mm <- turi2_df_library_DMSO_scalc_dcalc$m/turi2_df_library_DMSO_scalc_dcalc$s_val
#need to determine cutoff
###############################################################
#write function for other datasets?
#need to input the spread and filtered version of bac/library dfs -> the name of the dataframe should NOT contain "spread"
#ssmd_test <- function(bac, lib){
#   #bac$diff <- bac$"8" -  bac$"1"
#   #turi2_df_bac_DMSO_avg <- mean(bac$diff)
#   #lib$diff <- lib$"8" -  lib$"1"
#   lib$diff_neg <- lib$diff - turi2_df_bac_DMSO_avg
#   lib_bywells <- group_by(lib, well) %>% summarize(m = mean(diff_neg))
#   lib_bywells_test <- left_join(lib, lib_bywells, by = "well")    
#   lib_bywells_test$sub <- lib_bywells_test$diff_neg - lib_bywells_test$m
#   lib_bywells_test$sub2 <- (lib_bywells_test$sub)^2
#   lib_scalc <- group_by(lib_bywells_test, well) %>% summarize(sum = sum(sub2))
#   lib_scalc$s_val <- sqrt((lib_scalc$sum)/2)
#   lib_scalc_dcalc <- left_join(lib_scalc, lib_bywells_test, by = "well")    
#   lib_scalc_dcalc$mm <- lib_scalc_dcalc$m/lib_scalc_dcalc$s_val
#   lib_scalc_dcalc <- lib_scalc_dcalc[!duplicated(lib_scalc_dcalc$well),]
#   return(lib_scalc_dcalc)
# }

#turns out that using the filtered version for turi does make a difference on the results, since there are some wells without growth
#turi2_MeOH_ssmd_2 <- ssmd_test(turi2_df_bac_MeOH_spread, turi2_df_library_MeOH_spread)
# turi2_MeOH_ssmd <- ssmd_test(turi2_df_bac_MeOH, turi2_df_library_MeOH)
# turi2_H2O_ssmd <- ssmd_test(turi2_df_bac_H2O, turi2_df_library_H2O)
# turi2_DMSO_ssmd <- ssmd_test(turi2_df_bac_DMSO, turi2_df_library_DMSO)

#join with platemap?
turi2_H2O_ssmd_map <- left_join(turi2_H2O_ssmd, platemap_H2O_merge_form, by = "well")
#we need the chemical names, so use the manually edited file that essentially maps the messed up names with the regex names
# THERE HAS TO BE A BETTER WAY OF DOING THIS, LIKE ASSIGNING A CODE OR SMTH, BUT IM DUMB and have no foresight haha
turi_jw_key <- read.csv("turi_out_fuzzy_jw_edited_DO_NOT_EDIT.csv")
colnames(turi_jw_key) <- c("X", "solv.x", "well", "category", "solv", "dist")
turi2_H2O_ssmd_map <- left_join(turi2_H2O_ssmd_map, turi_jw_key, by = "solv") 
#assign color based on mm value
colfunc <- colorRampPalette(c('purple', 'green'))

#essentially using the order of the the "mm" value, arrange and assign a rank column with the row number
#testing for function########
turi2_H2O_ssmd_map_color <- turi2_H2O_ssmd_map %>%
  drop_na(mm) %>% 
  arrange(mm) %>%
  mutate(rank = row_number())

#join the initial dataframe and the dataframe with the color gradient by the mm value, and "rank" the colors generated 
#i'm not 100% sure how this works either
turi2_H2O_ssmd_map <- left_join(turi2_H2O_ssmd_map, turi2_H2O_ssmd_map_color, by = "mm")
the_colors <- tibble(rank = 1:max(turi2_H2O_ssmd_map$rank, na.rm = T), color_code = colfunc(max(turi2_H2O_ssmd_map$rank, na.rm = T)))

turi2_H2O_ssmd_map_color_test <- turi2_H2O_ssmd_map %>%
  left_join(the_colors, by = 'rank') %>%
  mutate(color_code = if_else(is.na(mm), '#FFFFFFFF', color_code))

#select out the columns of interest from the mess of columns that is a result of joining over and over again
turi2_H2O_ssmd_map_color_test <- turi2_H2O_ssmd_map_color_test[!duplicated(turi2_H2O_ssmd_map_color_test$solv.x.y),]
turi2_H2O_ssmd_map_color_test_clean <- select(turi2_H2O_ssmd_map_color_test, solv.x.y, color_code)

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
turi2_DMSO_color_out <- color_function(turi2_DMSO_ssmd, platemap_DMSO_merge)
turi2_MeOH_color_out <- color_function(turi2_MeOH_ssmd, platemap_MeOH_merge)
turi2_H2O_color_out <- color_function(turi2_H2O_ssmd, platemap_H2O_merge_form)
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
turi2_DMSO_ssmd_color_2 <- turi2_DMSO_ssmd
turi2_MeOH_ssmd_color_2 <- turi2_MeOH_ssmd
turi2_H2O_ssmd_color_2 <- turi2_H2O_ssmd


#using nested ifelse, find values above and below the cutoffs\
#there is probably a better way to do this using "when"
####note to self I think there is an extra step in here
turi2_DMSO_ssmd_color_2$color_2 <- ifelse(turi2_DMSO_ssmd_color_2$mm <= -1.645, "#8b12b3", 
                                          ifelse(turi2_DMSO_ssmd_color_2$mm >=1.645, "#12b31d",
                                                 ifelse(turi2_DMSO_ssmd_color_2$mm > -1.645 & turi2_DMSO_ssmd_color_2$mm < 1.645, "#8f8f8f", "NA")))
turi2_H2O_ssmd_color_2$color_2 <- ifelse(turi2_H2O_ssmd_color_2$mm <= -1.645, "#8b12b3", 
                                         ifelse(turi2_H2O_ssmd_color_2$mm >=1.645, "#12b31d",
                                                ifelse(turi2_H2O_ssmd_color_2$mm > -1.645 & turi2_H2O_ssmd_color_2$mm < 1.645, "#8f8f8f", "NA")))
turi2_MeOH_ssmd_color_2$color_2 <- ifelse(turi2_MeOH_ssmd_color_2$mm <= -1.645, "#8b12b3", 
                                          ifelse(turi2_MeOH_ssmd_color_2$mm >=1.645, "#12b31d",
                                                 ifelse(turi2_MeOH_ssmd_color_2$mm > -1.645 & turi2_MeOH_ssmd_color_2$mm < 1.645, "#8f8f8f", "NA")))


#correlate to chemical names from regex
turi2_DMSO_ssmd_color_2 <- left_join(turi2_DMSO_ssmd_color_2, platemap_DMSO_merge, by = "well")
turi2_DMSO_ssmd_color_2 <- select(turi2_DMSO_ssmd_color_2, well, mm, color_2, solv)
colnames(turi2_DMSO_ssmd_color_2) <- c("well", "mm", "color_2", "solv")

turi2_MeOH_ssmd_color_2 <- left_join(turi2_MeOH_ssmd_color_2, platemap_MeOH_merge, by = "well")
turi2_MeOH_ssmd_color_2 <- select(turi2_MeOH_ssmd_color_2, well, mm, color_2, solv)

turi2_H2O_ssmd_color_2 <- left_join(turi2_H2O_ssmd_color_2, platemap_H2O_merge_form, by = "well")
turi2_H2O_ssmd_color_2 <- select(turi2_H2O_ssmd_color_2, well, mm, color_2, solv)

#combine
test_turi2_color_2 <- rbind(turi2_DMSO_ssmd_color_2, turi2_H2O_ssmd_color_2, turi2_MeOH_ssmd_color_2)
table(test_turi2_color_2$color_2)
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


######try ssmd for splines maybe?

#calculate splines using lapply, put back into dataframe
#do ssmd method above, but on "spline" column isntead of "diff" column