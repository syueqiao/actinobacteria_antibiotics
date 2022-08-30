#script neatly combining entire methods pipeline for turicella

#load packages that are needed
library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(lspline)
library(fuzzyjoin)

#set wd
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/strains_data")
#set this function
`%notin%` <- Negate(`%in%`)

#import platemap
#see supporting files for this document
platemap <- read.csv("../plate_maps.csv", header = TRUE)

#################################TO DO

#####MAKE DMSO PLATEMAP FILES########################################
platemap_DMSO <- select(platemap, Well, DMSO)
DMSO_media_list <- platemap_DMSO[platemap_DMSO$DMSO %like% "Empty",]  
DMSO_media_list <- as.vector(DMSO_media_list$Well)

DMSO_bac_control_list <- platemap_DMSO[platemap_DMSO$DMSO %like% "Bacteria only", ]  
DMSO_bac_control_list <- as.vector(DMSO_bac_control_list$Well)

DMSO_positive_control_list <- platemap_DMSO[platemap_DMSO$DMSO %like% "Positive", ]
DMSO_positive_control_list <- as.vector(DMSO_positive_control_list$Well)
#create list of controls and stuff to filter out
DMSO_filter_out <- c(DMSO_media_list, DMSO_positive_control_list, DMSO_bac_control_list)

####MAKE H2O PLATEMAP FILES#########################################
platemap_H2O <- select(platemap, Well, H2O)
H2O_media_list <- platemap_H2O[platemap_H2O$H2O %like% "Empty",]  
H2O_media_list <- as.vector(H2O_media_list$Well)

H2O_bac_control_list <- platemap_H2O[platemap_H2O$H2O %like% "Bacteria only", ]  
H2O_bac_control_list <- as.vector(H2O_bac_control_list$Well)

H2O_positive_control_list <- platemap_H2O[platemap_H2O$H2O %like% "Positive", ]
H2O_positive_control_list <- as.vector(H2O_positive_control_list$Well)
#create list of controls and stuff to filter out
H2O_filter_out <- c(H2O_media_list, H2O_positive_control_list, H2O_bac_control_list)

#####MAKE MEOH PLATEMAP FILES####
platemap_MeOH <- select(platemap, Well, MeOH)
MeOH_media_list <- platemap_MeOH[platemap_MeOH$MeOH %like% "Empty",]  
MeOH_media_list <- as.vector(MeOH_media_list$Well)

MeOH_bac_control_list <- platemap_MeOH[platemap_MeOH$MeOH %like% "Bacteria only", ]  
MeOH_bac_control_list <- as.vector(MeOH_bac_control_list$Well)

MeOH_positive_control_list <- platemap_MeOH[platemap_MeOH$MeOH %like% "Positive", ]
MeOH_positive_control_list <- as.vector(MeOH_positive_control_list$Well)
#create list of controls and stuff to filter out
MeOH_filter_out <- c(MeOH_media_list, MeOH_positive_control_list, MeOH_bac_control_list)

###############set up loop for importing the metadata based ON THE DIRECTORY NAMES!!!
strains<- list.files(path="./", full.names = T , recursive =F)
solvents <- list.files(path=paste0(strains), full.names = T , recursive =F)
reps <- list.files(path=paste0(solvents), full.names = T, recursive = F)
files <- list.files(path=paste0(reps), pattern="*.CSV", full.names=TRUE, recursive=FALSE)
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


#create function to graph the curves in spread dataframes

all_curves_solvent_new <- function(df){
  col_num <- ncol(df)
  df <- select(df, well, rep, c(7:col_num))
  df$id <- paste(df$well, df$rep)
  col_num2 <- ncol(df)
  df <- select(df, c(3:col_num2))
  df <- melt(df, id.vars = 'id')
  
  
  return(ggplot(df, aes(x = variable, y = value)) + 
           geom_line(aes(color = id, group = id)) + 
           theme(legend.position="none"))
}

remove_first_time <- function(input){
  input_1 <- filter(input, time != "1")
  input_1$time <- as.numeric(input_1$time)
  input_1$time <- input_1$time-1 
  return(input_1)
}
#filter data for turicella only

#make inputs for turi2#####
turi2_df <- filter(df_metadata, bug == "turicella2")
#assign unique ID for each row (well and solvent)
turi2_df$UID <- paste(turi2_df$well, turi2_df$solvent, sep = "")
#separate into solvents for QC and stuff
turi2_df_DMSO <- filter(turi2_df , solvent == "DMSO")
turi2_df_H2O <- filter(turi2_df , solvent == "H2O")
turi2_df_MeOH <- filter(turi2_df , solvent == "MeOH")

#remove first two time points to account for the shaking time point and bubbles
turi2_df_library_DMSO <- filter(turi2_df_DMSO, well %notin% DMSO_filter_out)
turi2_df_library_DMSO <- remove_first_time(turi2_df_library_DMSO)
turi2_df_library_DMSO <- remove_first_time(turi2_df_library_DMSO)
turi2_df_library_DMSO_spread <- spread(turi2_df_library_DMSO, key = time, value = OD)
# all_curves_solvent_new(turi2_df_library_DMSO_spread)

#make dataframe with only the control wells
turi2_df_bac_DMSO <- filter(turi2_df_DMSO, well %in% DMSO_bac_control_list)
turi2_df_bac_DMSO <- remove_first_time(turi2_df_bac_DMSO)
turi2_df_bac_DMSO <- remove_first_time(turi2_df_bac_DMSO)
turi2_df_bac_DMSO_spread <- spread(turi2_df_bac_DMSO, key = time, value = OD)
# all_curves_solvent_new(turi2_df_bac_DMSO_spread)

#make dataframe with only the antibiot ic wells
turi2_df_pos_DMSO <- filter(turi2_df_DMSO, well %in% DMSO_positive_control_list)
turi2_df_pos_DMSO <- remove_first_time(turi2_df_pos_DMSO)
turi2_df_pos_DMSO <- remove_first_time(turi2_df_pos_DMSO)
turi2_df_pos_DMSO_spread <- spread(turi2_df_pos_DMSO, key = time, value = OD)
# all_curves_solvent_new(turi2_df_pos_DMSO_spread)

#make similar dataframe but for MeOH
turi2_df_library_MeOH <- filter(turi2_df_MeOH, well %notin% MeOH_filter_out)
turi2_df_library_MeOH <- remove_first_time(turi2_df_library_MeOH)
turi2_df_library_MeOH <- remove_first_time(turi2_df_library_MeOH)
turi2_df_library_MeOH_spread <- spread(turi2_df_library_MeOH, key = time, value = OD)
# all_curves_solvent_new(turi2_df_library_MeOH_spread)

#make dataframe with only the control wells
#only one with filtering, as there are some wells where thebacteria did not grow
turi2_df_bac_MeOH <- filter(turi2_df_MeOH, well %in% MeOH_bac_control_list)
turi2_df_bac_MeOH <- remove_first_time(turi2_df_bac_MeOH)
turi2_df_bac_MeOH <- remove_first_time(turi2_df_bac_MeOH)
turi2_df_bac_MeOH_spread <- spread(turi2_df_bac_MeOH, key = time, value = OD)
turi2_df_bac_MeOH_spread <- filter(turi2_df_bac_MeOH_spread, turi2_df_bac_MeOH_spread$"8" > 0.5)
# all_curves_solvent_new(turi2_df_bac_MeOH_spread)

#make dataframe with only the antibiotic wells
turi2_df_pos_MeOH <- filter(turi2_df_MeOH, well %in% MeOH_positive_control_list)
turi2_df_pos_MeOH <- remove_first_time(turi2_df_pos_MeOH)
turi2_df_pos_MeOH <- remove_first_time(turi2_df_pos_MeOH)
turi2_df_pos_MeOH_spread <- spread(turi2_df_pos_MeOH, key = time, value = OD)
# all_curves_solvent_new(turi2_df_pos_MeOH_spread)

#make for H2O
turi2_df_library_H2O <- filter(turi2_df_H2O, well %notin% H2O_filter_out)
turi2_df_library_H2O <- remove_first_time(turi2_df_library_H2O)
turi2_df_library_H2O <- remove_first_time(turi2_df_library_H2O)
turi2_df_library_H2O_spread <- spread(turi2_df_library_H2O, key = time, value = OD)
# all_curves_solvent_new(turi2_df_library_H2O_spread)

#make dataframe with only the control wells
turi2_df_bac_H2O <- filter(turi2_df_H2O, well %in% H2O_bac_control_list)
turi2_df_bac_H2O <- remove_first_time(turi2_df_bac_H2O)
turi2_df_bac_H2O <- remove_first_time(turi2_df_bac_H2O)
turi2_df_bac_H2O_spread <- spread(turi2_df_bac_H2O, key = time, value = OD)
# all_curves_solvent_new(turi2_df_bac_H2O_spread)

#make dataframe with only the antibiot ic wells
turi2_df_pos_H2O <- filter(turi2_df_H2O, well %in% H2O_positive_control_list)
turi2_df_pos_H2O <- remove_first_time(turi2_df_pos_H2O)
turi2_df_pos_H2O <- remove_first_time(turi2_df_pos_H2O)
turi2_df_pos_H2O_spread <- spread(turi2_df_pos_H2O, key = time, value = OD)
# all_curves_solvent_new(turi2_df_pos_H2O_spread)

#test the quality of each plate using SSMD
#define functions for calculating SSMDQC for spline and od600
#od600
ssmd_QC_diff <- function(df_bac, df_pos, df_rep){
  strain_df_bac_DMSO_rep <- filter(df_bac, rep == df_rep)
  #calc od600 and related values for bac and antibiotic control
  strain_df_bac_DMSO_rep$diff <- strain_df_bac_DMSO_rep$"8" - strain_df_bac_DMSO_rep$"1"
  strain_df_bac_DMSO_rep_mean <- mean(strain_df_bac_DMSO_rep$diff)
  strain_df_bac_DMSO_rep_var <- var(strain_df_bac_DMSO_rep$diff)
  
  strain_df_pos_DMSO_rep <- filter(df_pos, rep == df_rep)
  strain_df_pos_DMSO_rep$diff <- strain_df_pos_DMSO_rep$"8" - strain_df_pos_DMSO_rep$"1"
  strain_df_pos_DMSO_rep_mean <- mean(strain_df_pos_DMSO_rep$diff)
  strain_df_pos_DMSO_rep_var <- var(strain_df_pos_DMSO_rep$diff)
  
  strain_df_pos_DMSO_rep_ssmd <- (strain_df_pos_DMSO_rep_mean - strain_df_bac_DMSO_rep_mean)/sqrt((strain_df_pos_DMSO_rep_var) + (strain_df_bac_DMSO_rep_var))
  return(strain_df_pos_DMSO_rep_ssmd)
}

#spline 
ssmd_QC_spline <- function(df_bac, df_pos, df_rep_spline){
  
  strain_df_bac_DMSO_rep_spline <- filter(df_bac, rep == df_rep_spline)
  
  #a bit more complicated, calculate the spline on the spot
  knots <- c(3)
  strain_df_bac_DMSO_rep_spline_values <- select(strain_df_bac_DMSO_rep_spline, c(7:14))
  time <- c(0:7)
  strain_df_bac_DMSO_rep_spline_vec <- apply(strain_df_bac_DMSO_rep_spline_values, 1, function(x) lm(x ~ lspline(time, knots = knots)))
  strain_df_bac_DMSO_rep_spline_vec_coeff <- as.data.frame(unlist(lapply(strain_df_bac_DMSO_rep_spline_vec, function(x) coef(x)[2])))
  strain_df_bac_DMSO_rep_spline$spline <- strain_df_bac_DMSO_rep_spline_vec_coeff
  
  strain_df_bac_DMSO_rep_spline_mean <- mean(strain_df_bac_DMSO_rep_spline$spline)
  strain_df_bac_DMSO_rep_spline_var <- var(strain_df_bac_DMSO_rep_spline$spline)
  
  #for postive control
  strain_df_pos_DMSO_rep_spline <- filter(df_pos, rep == df_rep_spline)
  
  strain_df_pos_DMSO_rep_spline_values <- select(strain_df_pos_DMSO_rep_spline, c(7:14))
  time <- c(0:7)
  strain_df_pos_DMSO_rep_spline_vec <- apply(strain_df_pos_DMSO_rep_spline_values, 1, function(x) lm(x ~ lspline(time, knots = knots)))
  strain_df_pos_DMSO_rep_spline_vec_coeff <- as.data.frame(unlist(lapply(strain_df_pos_DMSO_rep_spline_vec, function(x) coef(x)[2])))
  strain_df_pos_DMSO_rep_spline$spline <- strain_df_pos_DMSO_rep_spline_vec_coeff
  
  strain_df_pos_DMSO_rep_spline_mean <- mean(strain_df_pos_DMSO_rep_spline$spline)
  strain_df_pos_DMSO_rep_spline_var <- var(strain_df_pos_DMSO_rep_spline$spline)
  
  strain_df_pos_DMSO_rep_spline_ssmd <- (strain_df_pos_DMSO_rep_spline_mean - strain_df_bac_DMSO_rep_spline_mean)/sqrt((strain_df_pos_DMSO_rep_spline_var) + (strain_df_bac_DMSO_rep_spline_var))
  return(strain_df_pos_DMSO_rep_spline_ssmd)}

#run the functions, record output values
#for the spline
ssmd_QC_spline(turi2_df_bac_H2O_spread, turi2_df_pos_H2O_spread, "rep1")
ssmd_QC_spline(turi2_df_bac_H2O_spread, turi2_df_pos_H2O_spread, "rep2")
ssmd_QC_spline(turi2_df_bac_H2O_spread, turi2_df_pos_H2O_spread, "rep3")

ssmd_QC_spline(turi2_df_bac_MeOH_spread, turi2_df_pos_MeOH_spread, "rep1")
ssmd_QC_spline(turi2_df_bac_MeOH_spread, turi2_df_pos_MeOH_spread, "rep2")
ssmd_QC_spline(turi2_df_bac_MeOH_spread, turi2_df_pos_MeOH_spread, "rep3")

ssmd_QC_spline(turi2_df_bac_DMSO_spread, turi2_df_pos_DMSO_spread, "rep1")
ssmd_QC_spline(turi2_df_bac_DMSO_spread, turi2_df_pos_DMSO_spread, "rep2")
ssmd_QC_spline(turi2_df_bac_DMSO_spread, turi2_df_pos_DMSO_spread, "rep3")

#for od600
ssmd_QC_diff(turi2_df_bac_H2O_spread, turi2_df_pos_H2O_spread, "rep1")
ssmd_QC_diff(turi2_df_bac_H2O_spread, turi2_df_pos_H2O_spread, "rep2")
ssmd_QC_diff(turi2_df_bac_H2O_spread, turi2_df_pos_H2O_spread, "rep3")

ssmd_QC_diff(turi2_df_bac_MeOH_spread, turi2_df_pos_MeOH_spread, "rep1")
ssmd_QC_diff(turi2_df_bac_MeOH_spread, turi2_df_pos_MeOH_spread, "rep2")
ssmd_QC_diff(turi2_df_bac_MeOH_spread, turi2_df_pos_MeOH_spread, "rep3")

ssmd_QC_diff(turi2_df_bac_DMSO_spread, turi2_df_pos_DMSO_spread, "rep1")
ssmd_QC_diff(turi2_df_bac_DMSO_spread, turi2_df_pos_DMSO_spread, "rep2")
ssmd_QC_diff(turi2_df_bac_DMSO_spread, turi2_df_pos_DMSO_spread, "rep3")

#next do the analysis 
#do everything for the SSMD method first, all the way up to getting categories

#troubleshooting -> uncomment below
# bac_in <- turi2_df_bac_DMSO_spread
# lib_in <- turi2_df_library_DMSO_spread

#define function that calculates the ssmd value using the method-of-moment method

OD600_ssmd <- function(bac_in, lib_in){
  #calculate dOD600
  bac_in$diff <- bac_in$"8" -  bac_in$"1"
  #separate into different reps to do the paired ssmd test
  turi2_df_bac_in_solv_avg_rep1 <- bac_in %>%
    filter(rep == "rep1") %>%
    .$diff %>%
    mean
  
  turi2_df_bac_in_solv_avg_rep2 <- bac_in %>%
    filter(rep == "rep2") %>%
    .$diff %>%
    mean
  
  turi2_df_bac_in_solv_avg_rep3 <- bac_in %>%
    filter(rep == "rep3") %>%
    .$diff %>%
    mean
  
  #repeat same process for the library plates
  lib_in$diff <- lib_in$"8" -  lib_in$"1"
  
  #start calculation for di, for each well
  turi2_df_lib_in_solv_avg_rep1 <- lib_in %>%
    filter(rep == "rep1")
  turi2_df_lib_in_solv_avg_rep1$d_1 <- turi2_df_lib_in_solv_avg_rep1$diff -  turi2_df_bac_in_solv_avg_rep1
  
  turi2_df_lib_in_solv_avg_rep2 <- lib_in %>%
    filter(rep == "rep2")
  turi2_df_lib_in_solv_avg_rep2$d_1 <- turi2_df_lib_in_solv_avg_rep2$diff -  turi2_df_bac_in_solv_avg_rep2
  
  turi2_df_lib_in_solv_avg_rep3 <- lib_in %>%
    filter(rep == "rep3")
  turi2_df_lib_in_solv_avg_rep3$d_1 <- turi2_df_lib_in_solv_avg_rep3$diff -  turi2_df_bac_in_solv_avg_rep3
  
  #combine into one df so you can average the values
  turi2_df_lib_in_solv_avg_repall <- rbind(turi2_df_lib_in_solv_avg_rep1, turi2_df_lib_in_solv_avg_rep2, turi2_df_lib_in_solv_avg_rep3)
 
  #group by wells, average those values
   turi2_df_lib_in_solv_avg_repall_avg <- group_by(turi2_df_lib_in_solv_avg_repall, well) %>%
    summarize(d_calc = mean(d_1)) %>%
    left_join(., turi2_df_lib_in_solv_avg_repall, by = "well")
  
  #calculate s -> find the variance of the d_1
  turi2_df_lib_in_solv_avg_repall_ds <- group_by(turi2_df_lib_in_solv_avg_repall_avg, well) %>%
    summarize(s_calc = sd(d_1)) %>%
    left_join(., turi2_df_lib_in_solv_avg_repall_avg, by = "well")
  
  #finally calculate ssmd
  turi2_df_lib_in_solv_avg_repall_ds$mm_diff <- turi2_df_lib_in_solv_avg_repall_ds$d_calc/turi2_df_lib_in_solv_avg_repall_ds$s_calc
  turi2_df_lib_in_solv_avg_repall_ds <- turi2_df_lib_in_solv_avg_repall_ds[!duplicated(turi2_df_lib_in_solv_avg_repall_ds$well),]
  return(turi2_df_lib_in_solv_avg_repall_ds)
}

#run function for each solvent plate separately
turi2_DMSO_ssmd_OD600 <- OD600_ssmd(turi2_df_bac_DMSO_spread, turi2_df_library_DMSO_spread)
turi2_MeOH_ssmd_OD600 <- OD600_ssmd(turi2_df_bac_MeOH_spread, turi2_df_library_MeOH_spread)
turi2_H2O_ssmd_OD600 <- OD600_ssmd(turi2_df_bac_H2O_spread, turi2_df_library_H2O_spread)

#assign unique ID so naming isn't so painful later
turi2_all_ssmd_OD600 <- rbind(turi2_DMSO_ssmd_OD600, turi2_MeOH_ssmd_OD600, turi2_H2O_ssmd_OD600)
turi2_all_ssmd_OD600$UID <- paste(turi2_all_ssmd_OD600$solvent, turi2_all_ssmd_OD600$well, sep = "")


#next, define function for doing the same ssmd test on the splines this time
# #troubleshooting
# bac_in_spline <- turi2_df_bac_DMSO_spread
# lib_in_spline <- turi2_df_library_DMSO_spread

spline_ssmd <- function(bac_in_spline, lib_in_spline){
  
  #subset values, for lspline
  turi2_df_bac_solv_splines <- select(bac_in_spline, c(7:14))
  #set vector for time
  turi2_bac_time_solv <- c(0:7)
  #set knots
  knots <- c(3)
  #make list with lspline
  turi2_bac_lm_solv <- apply(turi2_df_bac_solv_splines, 1, function(x) lm(x ~ lspline(turi2_bac_time_solv, knots = knots)))
  #choose the spline
  turi2_df_bac_solv_splines_coeff <- as.data.frame(unlist(lapply(turi2_bac_lm_solv, function(x) coef(x)[2])))
  
  #add spline value to df , follow same workflow as OD600, since we are just applying the same stuff to different values
  bac_in_spline$spline <- turi2_df_bac_solv_splines_coeff

  turi2_df_bac_in_spline_solv_avg_rep1 <- bac_in_spline %>%
    filter(rep == "rep1") %>%
    .$spline %>%
    mean
  
  turi2_df_bac_in_spline_solv_avg_rep2 <- bac_in_spline %>%
    filter(rep == "rep2") %>%
    .$spline %>%
    mean
  
  turi2_df_bac_in_spline_solv_avg_rep3 <- bac_in_spline %>%
    filter(rep == "rep3") %>%
    .$spline %>%
    mean
  
  #calculate spline for library wells now, following same logic as above
  turi2_df_library_solv_splines <- select(lib_in_spline, c(7:14))
  turi2_library_time_solv <- c(0:7)
  turi2_library_lm_solv <- apply(turi2_df_library_solv_splines, 1, function(x) lm(x ~ lspline(turi2_library_time_solv, knots = knots)))
  turi2_df_library_solv_splines_coeff <- as.data.frame(unlist(lapply(turi2_library_lm_solv, function(x) coef(x)[2])))
  

  lib_in_spline$spline <- turi2_df_library_solv_splines_coeff

  turi2_df_lib_in_spline_solv_avg_rep1 <- lib_in_spline %>%
    filter(rep == "rep1")
  turi2_df_lib_in_spline_solv_avg_rep1$d_1 <- turi2_df_lib_in_spline_solv_avg_rep1$spline -  turi2_df_bac_in_spline_solv_avg_rep1
  
  turi2_df_lib_in_spline_solv_avg_rep2 <- lib_in_spline %>%
    filter(rep == "rep2")
  turi2_df_lib_in_spline_solv_avg_rep2$d_1 <- turi2_df_lib_in_spline_solv_avg_rep2$spline -  turi2_df_bac_in_spline_solv_avg_rep2
  
  turi2_df_lib_in_spline_solv_avg_rep3 <- lib_in_spline %>%
    filter(rep == "rep3")
  turi2_df_lib_in_spline_solv_avg_rep3$d_1 <- turi2_df_lib_in_spline_solv_avg_rep3$spline -  turi2_df_bac_in_spline_solv_avg_rep3
  
  turi2_df_lib_in_spline_solv_avg_repall <- rbind(turi2_df_lib_in_spline_solv_avg_rep1, turi2_df_lib_in_spline_solv_avg_rep2, turi2_df_lib_in_spline_solv_avg_rep3)
  
  turi2_df_lib_in_spline_solv_avg_repall_avg <- group_by(turi2_df_lib_in_spline_solv_avg_repall, well) %>%
    summarize(d_calc = mean(d_1)) %>%
    left_join(., turi2_df_lib_in_spline_solv_avg_repall, by = "well")
  
  turi2_df_lib_in_spline_solv_avg_repall_ds <- group_by(turi2_df_lib_in_spline_solv_avg_repall_avg, well) %>%
    summarize(s_calc = sd(d_1)) %>%
    left_join(., turi2_df_lib_in_spline_solv_avg_repall_avg, by = "well")
  
  turi2_df_lib_in_spline_solv_avg_repall_ds$mm_spline <- turi2_df_lib_in_spline_solv_avg_repall_ds$d_calc/turi2_df_lib_in_spline_solv_avg_repall_ds$s_calc
  turi2_df_lib_in_spline_solv_avg_repall_ds <- turi2_df_lib_in_spline_solv_avg_repall_ds[!duplicated(turi2_df_lib_in_spline_solv_avg_repall_ds$well),]
  return(turi2_df_lib_in_spline_solv_avg_repall_ds)
}

turi2_DMSO_ssmd_spline <- spline_ssmd(turi2_df_bac_DMSO_spread, turi2_df_library_DMSO_spread)
turi2_MeOH_ssmd_spline <- spline_ssmd(turi2_df_bac_MeOH_spread, turi2_df_library_MeOH_spread)
turi2_H2O_ssmd_spline <- spline_ssmd(turi2_df_bac_H2O_spread, turi2_df_library_H2O_spread)

turi2_all_ssmd_spline <- rbind(turi2_DMSO_ssmd_spline, turi2_MeOH_ssmd_spline, turi2_H2O_ssmd_spline)
turi2_all_ssmd_spline$UID <- paste(turi2_all_ssmd_spline$solvent, turi2_all_ssmd_spline$well, sep = "")

#lean up before merging
turi2_all_ssmd_spline_trim <- select(turi2_all_ssmd_spline, UID, mm_spline)

turi2_all_ssmd_both <- left_join(turi2_all_ssmd_spline_trim, turi2_all_ssmd_OD600, by = "UID")
#assign categories
#first for od600
turi2_all_ssmd_both$cat_diff <- ifelse(turi2_all_ssmd_both$mm_diff <= -1.645, "inhib", 
                                            ifelse(turi2_all_ssmd_both$mm_diff >=1.645, "enh",
                                                   ifelse(turi2_all_ssmd_both$mm_diff > -1.645 & turi2_all_ssmd_both$mm_diff < 1.645, "noef", "NA")))
#next for spline
turi2_all_ssmd_both$cat_spline <- ifelse(turi2_all_ssmd_both$mm_spline <= -1.645, "inhib", 
                                         ifelse(turi2_all_ssmd_both$mm_spline >=1.645, "enh",
                                                ifelse(turi2_all_ssmd_both$mm_spline > -1.645 & turi2_all_ssmd_both$mm_spline < 1.645, "noef", "NA")))

#define function for categorizing

category_func_ssmd <- function(df){
  df$ssmd_category <- ifelse(df$cat_diff == "noef" & df$cat_spline == "noef", 5, 
                        ifelse(df$cat_diff == "inhib" & df$cat_spline == "inhib", 9, 
                               ifelse(df$cat_diff == "enh" & df$cat_spline == "enh", 1, 
                                      ifelse(df$cat_diff == "enh" & df$cat_spline == "noef", 2, 
                                             ifelse(df$cat_diff == "enh" & df$cat_spline == "inhib", 3, 
                                                    ifelse(df$cat_diff == "noef" & df$cat_spline == "enh", 4,
                                                           ifelse(df$cat_diff == "noef" & df$cat_spline == "inhib", 6, 
                                                                  ifelse(df$cat_diff == "inhib" & df$cat_spline == "enh", 7,
                                                                         ifelse(df$cat_diff == "inhib" & df$cat_spline == "noef", 8, 0)))))))))
  hist(df$ssmd_category)
  return(df)
}

turi2_all_ssmd_both_category <- category_func_ssmd(turi2_all_ssmd_both)

#do od600 and spline for the t test now


OD600_ttest <- function(bac_in, lib_in){
  #calculate dOD600
  bac_in$diff <- bac_in$"8" -  bac_in$"1"
  #separate into different reps to do the paired ttest test
  turi2_df_bac_in_solv_avg_rep1 <- bac_in %>%
    filter(rep == "rep1") %>%
    .$diff %>%
    mean
  
  turi2_df_bac_in_solv_avg_rep2 <- bac_in %>%
    filter(rep == "rep2") %>%
    .$diff %>%
    mean
  
  turi2_df_bac_in_solv_avg_rep3 <- bac_in %>%
    filter(rep == "rep3") %>%
    .$diff %>%
    mean
  
  #repeat same process for the library plates
  lib_in$diff <- lib_in$"8" -  lib_in$"1"
  
  #start calculation for di, for each well
  turi2_df_lib_in_solv_avg_rep1 <- lib_in %>%
    filter(rep == "rep1")
  turi2_df_lib_in_solv_avg_rep1$d_1 <- turi2_df_lib_in_solv_avg_rep1$diff -  turi2_df_bac_in_solv_avg_rep1
  
  turi2_df_lib_in_solv_avg_rep2 <- lib_in %>%
    filter(rep == "rep2")
  turi2_df_lib_in_solv_avg_rep2$d_1 <- turi2_df_lib_in_solv_avg_rep2$diff -  turi2_df_bac_in_solv_avg_rep2
  
  turi2_df_lib_in_solv_avg_rep3 <- lib_in %>%
    filter(rep == "rep3")
  turi2_df_lib_in_solv_avg_rep3$d_1 <- turi2_df_lib_in_solv_avg_rep3$diff -  turi2_df_bac_in_solv_avg_rep3
  
  #combine into one df so you can average the values
  turi2_df_lib_in_solv_avg_repall <- rbind(turi2_df_lib_in_solv_avg_rep1, turi2_df_lib_in_solv_avg_rep2, turi2_df_lib_in_solv_avg_rep3)
  
  #group by wells, average those values
  turi2_df_lib_in_solv_avg_repall_avg <- group_by(turi2_df_lib_in_solv_avg_repall, well) %>%
    summarize(d_calc = mean(d_1)) %>%
    left_join(., turi2_df_lib_in_solv_avg_repall, by = "well")
  
  #calculate s -> find the variance of the d_1
  turi2_df_lib_in_solv_avg_repall_ds <- group_by(turi2_df_lib_in_solv_avg_repall_avg, well) %>%
    summarize(s_calc = sd(d_1)) %>%
    left_join(., turi2_df_lib_in_solv_avg_repall_avg, by = "well")
  
  #finally calculate ttest
  turi2_df_lib_in_solv_avg_repall_ds$tt_diff <- (sqrt(3)*turi2_df_lib_in_solv_avg_repall_ds$d_calc)/turi2_df_lib_in_solv_avg_repall_ds$s_calc
  turi2_df_lib_in_solv_avg_repall_ds <- turi2_df_lib_in_solv_avg_repall_ds[!duplicated(turi2_df_lib_in_solv_avg_repall_ds$well),]
  turi2_df_lib_in_solv_avg_repall_ds$p_value <- 2*pt(q = abs(turi2_df_lib_in_solv_avg_repall_ds$tt_diff), df = 2, lower.tail = F)
  return(turi2_df_lib_in_solv_avg_repall_ds)
}

#run function for each solvent plate separately
turi2_DMSO_ttest_OD600 <- OD600_ttest(turi2_df_bac_DMSO_spread, turi2_df_library_DMSO_spread)
turi2_MeOH_ttest_OD600 <- OD600_ttest(turi2_df_bac_MeOH_spread, turi2_df_library_MeOH_spread)
turi2_H2O_ttest_OD600 <- OD600_ttest(turi2_df_bac_H2O_spread, turi2_df_library_H2O_spread)

#assign unique ID so naming isn't so painful later
turi2_all_ttest_OD600 <- rbind(turi2_DMSO_ttest_OD600, turi2_MeOH_ttest_OD600, turi2_H2O_ttest_OD600)
turi2_all_ttest_OD600$UID <- paste(turi2_all_ttest_OD600$solvent, turi2_all_ttest_OD600$well, sep = "")


#next, define function for doing the same ttest test on the splines this time
# #troubleshooting
# bac_in_spline <- turi2_df_bac_DMSO_spread
# lib_in_spline <- turi2_df_library_DMSO_spread

spline_ttest <- function(bac_in_spline, lib_in_spline){
  
  #subset values, for lspline
  turi2_df_bac_solv_splines <- select(bac_in_spline, c(7:14))
  #set vector for time
  turi2_bac_time_solv <- c(0:7)
  #set knots
  knots <- c(3)
  #make list with lspline
  turi2_bac_lm_solv <- apply(turi2_df_bac_solv_splines, 1, function(x) lm(x ~ lspline(turi2_bac_time_solv, knots = knots)))
  #choose the spline
  turi2_df_bac_solv_splines_coeff <- as.data.frame(unlist(lapply(turi2_bac_lm_solv, function(x) coef(x)[2])))
  
  #add spline value to df , follow same workflow as OD600, since we are just applying the same stuff to different values
  bac_in_spline$spline <- turi2_df_bac_solv_splines_coeff
  
  turi2_df_bac_in_spline_solv_avg_rep1 <- bac_in_spline %>%
    filter(rep == "rep1") %>%
    .$spline %>%
    mean
  
  turi2_df_bac_in_spline_solv_avg_rep2 <- bac_in_spline %>%
    filter(rep == "rep2") %>%
    .$spline %>%
    mean
  
  turi2_df_bac_in_spline_solv_avg_rep3 <- bac_in_spline %>%
    filter(rep == "rep3") %>%
    .$spline %>%
    mean
  
  #calculate spline for library wells now, following same logic as above
  turi2_df_library_solv_splines <- select(lib_in_spline, c(7:14))
  turi2_library_time_solv <- c(0:7)
  turi2_library_lm_solv <- apply(turi2_df_library_solv_splines, 1, function(x) lm(x ~ lspline(turi2_library_time_solv, knots = knots)))
  turi2_df_library_solv_splines_coeff <- as.data.frame(unlist(lapply(turi2_library_lm_solv, function(x) coef(x)[2])))
  
  
  lib_in_spline$spline <- turi2_df_library_solv_splines_coeff
  
  turi2_df_lib_in_spline_solv_avg_rep1 <- lib_in_spline %>%
    filter(rep == "rep1")
  turi2_df_lib_in_spline_solv_avg_rep1$d_1 <- turi2_df_lib_in_spline_solv_avg_rep1$spline -  turi2_df_bac_in_spline_solv_avg_rep1
  
  turi2_df_lib_in_spline_solv_avg_rep2 <- lib_in_spline %>%
    filter(rep == "rep2")
  turi2_df_lib_in_spline_solv_avg_rep2$d_1 <- turi2_df_lib_in_spline_solv_avg_rep2$spline -  turi2_df_bac_in_spline_solv_avg_rep2
  
  turi2_df_lib_in_spline_solv_avg_rep3 <- lib_in_spline %>%
    filter(rep == "rep3")
  turi2_df_lib_in_spline_solv_avg_rep3$d_1 <- turi2_df_lib_in_spline_solv_avg_rep3$spline -  turi2_df_bac_in_spline_solv_avg_rep3
  
  turi2_df_lib_in_spline_solv_avg_repall <- rbind(turi2_df_lib_in_spline_solv_avg_rep1, turi2_df_lib_in_spline_solv_avg_rep2, turi2_df_lib_in_spline_solv_avg_rep3)
  
  turi2_df_lib_in_spline_solv_avg_repall_avg <- group_by(turi2_df_lib_in_spline_solv_avg_repall, well) %>%
    summarize(d_calc = mean(d_1)) %>%
    left_join(., turi2_df_lib_in_spline_solv_avg_repall, by = "well")
  
  turi2_df_lib_in_spline_solv_avg_repall_ds <- group_by(turi2_df_lib_in_spline_solv_avg_repall_avg, well) %>%
    summarize(s_calc = sd(d_1)) %>%
    left_join(., turi2_df_lib_in_spline_solv_avg_repall_avg, by = "well")
  
  turi2_df_lib_in_spline_solv_avg_repall_ds$tt_spline <- (sqrt(3)*turi2_df_lib_in_spline_solv_avg_repall_ds$d_calc)/turi2_df_lib_in_spline_solv_avg_repall_ds$s_calc
  turi2_df_lib_in_spline_solv_avg_repall_ds <- turi2_df_lib_in_spline_solv_avg_repall_ds[!duplicated(turi2_df_lib_in_spline_solv_avg_repall_ds$well),]
  turi2_df_lib_in_spline_solv_avg_repall_ds$p_value_spline <- 2*pt(q = abs(turi2_df_lib_in_spline_solv_avg_repall_ds$tt_spline), df = 2, lower.tail = F)
  return(turi2_df_lib_in_spline_solv_avg_repall_ds)
}

turi2_DMSO_ttest_spline <- spline_ttest(turi2_df_bac_DMSO_spread, turi2_df_library_DMSO_spread)
turi2_MeOH_ttest_spline <- spline_ttest(turi2_df_bac_MeOH_spread, turi2_df_library_MeOH_spread)
turi2_H2O_ttest_spline <- spline_ttest(turi2_df_bac_H2O_spread, turi2_df_library_H2O_spread)

turi2_all_ttest_spline <- rbind(turi2_DMSO_ttest_spline, turi2_MeOH_ttest_spline, turi2_H2O_ttest_spline)
turi2_all_ttest_spline$UID <- paste(turi2_all_ttest_spline$solvent, turi2_all_ttest_spline$well, sep = "")

#lean up before merging
turi2_all_ttest_spline_trim <- select(turi2_all_ttest_spline, UID, tt_spline, p_value_spline)

turi2_all_ttest_both <- left_join(turi2_all_ttest_spline_trim, turi2_all_ttest_OD600, by = "UID")
#assign categories
#first for od600
turi2_all_ttest_both$cat_diff <- ifelse(turi2_all_ttest_both$tt_diff < 0 & turi2_all_ttest_both$p_value < 0.05, "inhib", 
                                       ifelse(turi2_all_ttest_both$p_value > 0.05, "noef",
                                              ifelse(turi2_all_ttest_both$tt_diff > 0 & turi2_all_ttest_both$p_value < 0.05, "enh", "NA")))
#next for spline
turi2_all_ttest_both$cat_spline <- ifelse(turi2_all_ttest_both$tt_spline < 0 & turi2_all_ttest_both$p_value_spline < 0.05, "inhib", 
                                          ifelse(turi2_all_ttest_both$p_value_spline > 0.05, "noef",
                                                 ifelse(turi2_all_ttest_both$tt_spline > 0 & turi2_all_ttest_both$p_value_spline < 0.05, "enh", "NA")))


category_func_ttest <- function(df){
  df$ttest_category <- ifelse(df$cat_diff == "noef" & df$cat_spline == "noef", 5, 
                             ifelse(df$cat_diff == "inhib" & df$cat_spline == "inhib", 9, 
                                    ifelse(df$cat_diff == "enh" & df$cat_spline == "enh", 1, 
                                           ifelse(df$cat_diff == "enh" & df$cat_spline == "noef", 2, 
                                                  ifelse(df$cat_diff == "enh" & df$cat_spline == "inhib", 3, 
                                                         ifelse(df$cat_diff == "noef" & df$cat_spline == "enh", 4,
                                                                ifelse(df$cat_diff == "noef" & df$cat_spline == "inhib", 6, 
                                                                       ifelse(df$cat_diff == "inhib" & df$cat_spline == "enh", 7,
                                                                              ifelse(df$cat_diff == "inhib" & df$cat_spline == "noef", 8, 0)))))))))
  hist(df$ttest_category)
  return(df)
}
turi2_all_ttest_both_category <- category_func_ttest(turi2_all_ttest_both)

#now, generate the chemical similarity tree for Sweet library
#to do this, need to do outside of R
#basically, take all the chemical names from the sweet library, then run it through the REGEX script that removes all the weird formatting
#run all the chemical names through pubchem identifier exchange service, changing from chemical names (synonyms) to SMILEs

#move output to txt file, add column headings, then import as tsvs
sweet_smiles <- read_tsv("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/chemicals_dupes.txt")

#remove duplicates because pubchem search returns duplicates (to be consistent, this function always keeps the first one)
sweet_smiles <- sweet_smiles[!duplicated(sweet_smiles$Chemical),]

#subset out the chemicals pubchem was not able to find smiles for
sweet_smiles_NA <- sweet_smiles[is.na(sweet_smiles$SMILE),]
write_tsv(sweet_smiles_NA, "sweet_smiles_NA.tsv")

#manually edit in the smiles that pubchem was not able to get
sweet_smiles_NA <- read_tsv("../sweet_smiles_NA.tsv")
sweet_smiles <- na.omit(sweet_smiles)

#combine two dataframes together 
sweet_smiles_all <- rbind(sweet_smiles, sweet_smiles_NA)
#reverse dataframes since chemmineR assumes the first column is smiles, and the second column is the identifier
sweet_smiles_all <- rev(sweet_smiles_all)
#write out to import
write_tsv(sweet_smiles_all, "sweet_smiles_import.txt")
#manually remove column names (blame ChemmineR)

#next, import the ChemmineR environment

smiset <- read.SMIset("sweet_smiles_import.txt")
#convert to chemminer working format
sdfset <- smiles2sdf(smiset)
#convert to apset format for clustering
apset <- sdf2ap(sdfset)
#important!!! run this to make sure that the number of compounds is what you expect
apset

#this is the clustering command from ChemmineR. I'm not exactly sure what the cutoffs do
#maybe take a look?
#this will save the file "distmat.rda" into your working directory
clusters <- cmp.cluster(db=apset, cutoff = c(0.7, 0.8, 0.9),
                        save.distances="distmat.rda")

#put distamat.rda into R environment
load("distmat.rda") 

#convert to matrix format
sweet_matrix <- as.matrix(distmat)

#ward.D2 uses non-squared distances in euclidean space. the distance matrix I have is not euclidean
#however, since we are not analyzing clusters rigourously and more for visualization purposes, should be okay
#previously erroneously used "ward" which should have squared values input
hc <- hclust(as.dist(distmat), method="ward.D2") 
hc[["labels"]] <- cid(apset) # Assign correct item labels 

#assign the hierarchal clustering output as a phylogenetic tree
#write out as a tree
t <- as.phylo(hc)
write.tree(t, "sweet_tree.txt")

#this output file then is uploaded to iTOL to visualize chemical structure similarities

#next step is to make the annotation file, to show which categories are which for each compound, for the ssmd and t test separately

#first, have to match the labels that change going into iTOL with the ones that we have on the UID map

#generating the UID map

platemap_DMSO_id <- platemap_DMSO
platemap_DMSO_id$UID <- paste("DMSO", platemap_DMSO$Well, sep = "")
platemap_DMSO_id <- filter(platemap_DMSO_id, Well %notin% DMSO_filter_out)
colnames(platemap_DMSO_id) <- c("well", "Chemical", "UID")

platemap_MeOH_id <- platemap_MeOH
platemap_MeOH_id$UID <- paste("MeOH", platemap_MeOH$Well, sep = "")
platemap_MeOH_id <- filter(platemap_MeOH_id, Well %notin% MeOH_filter_out)
colnames(platemap_MeOH_id) <- c("well", "Chemical", "UID")

platemap_H2O_id <- platemap_H2O
platemap_H2O_id$UID <- paste("H2O", platemap_H2O$Well, sep = "")
platemap_H2O_id <- filter(platemap_H2O_id, Well %notin% H2O_filter_out)
colnames(platemap_H2O_id) <- c("well", "Chemical", "UID")


platemap_id <- rbind(platemap_DMSO_id, platemap_MeOH_id, platemap_H2O_id)
platemap_id$well <- NULL

#combine regex names with OG names
platemap_id_key <- stringdist_join(platemap_id, sweet_smiles_all, method = "jw", by = "Chemical", max_dist = 99, distance_col = "dist") %>% 
  group_by(Chemical.x) %>% 
  slice_min(order_by = dist, n=1)

#check if columns have duplicates/if it worked
duplicated(platemap_id_key$Chemical.x)
#nice!
colnames(platemap_id_key) <- c("chemical_og", "UID", "SMILE", "chemical_regex", "dist")
#don't need this colum  anymore
platemap_id_key$dist <- NULL

#now we want to correlate the iTOL labels to the UIDs
#the iTOL format removes the first character if it is "(" -> we can accomplish this using simple regex
platemap_id_key$chemical_iTOL <- str_remove(platemap_id_key$chemical_regex, "^\\(+")
platemap_id_key$chemical_iTOL <- str_remove(platemap_id_key$chemical_iTOL, "\\)+$")
#itol format also changes "(" to "_"
platemap_id_key$chemical_iTOL <- gsub('\\(', '-',
                                      gsub('\\)', '-', 
                                           gsub(" ", "", 
                                                gsub(",", "-", 
                                                     gsub(":", "-", platemap_id_key$chemical_iTOL)))))


#select the relevant columns, the iTOL name, and the UID for matching to the categories
platemap_id_key_output <- select(ungroup(platemap_id_key), UID, chemical_iTOL)

#set new dfs so we don't overwrite them
turi2_ttest <- turi2_all_ttest_both_category
turi2_ssmd <- turi2_all_ssmd_both_category

#combine the iTOL chemical names with these df using the UID
turi2_ttest_annotation <- left_join(turi2_ttest, platemap_id_key_output, by = "UID")
turi2_ssmd_annotation <- left_join(turi2_ssmd, platemap_id_key_output, by = "UID")

#add colors to the dataframe based on the category that was assigned
turi2_ttest_annotation$ttest_color <- ifelse(turi2_ttest_annotation$ttest_category >= 1 & turi2_ttest_annotation$ttest_category <= 3, "#07eb29", 
                                          ifelse(turi2_ttest_annotation$ttest_category >= 4 & turi2_ttest_annotation$ttest_category <= 6 , "#a1a1a1",
                                                 ifelse(turi2_ttest_annotation$ttest_category >= 7 & turi2_ttest_annotation$ttest_category <= 9, "#ad05e6", "ur bad")))


turi2_ssmd_annotation$ssmd_color <- ifelse(turi2_ssmd_annotation$ssmd_category >= 1 & turi2_ssmd_annotation$ssmd_category <= 3, "#07eb29", 
                                             ifelse(turi2_ssmd_annotation$ssmd_category >= 4 & turi2_ssmd_annotation$ssmd_category <= 6 , "#a1a1a1",
                                                    ifelse(turi2_ssmd_annotation$ssmd_category >= 7 & turi2_ssmd_annotation$ssmd_category <= 9, "#ad05e6", "ur bad")))
#format the dataframe to the layout that iTOL needs
#separately for ssmd and ttest
turi2_ttest_annotation$col1 <- 2
turi2_ttest_annotation$col2 <- 2
turi2_ttest_annotation$col3 <- 1
turi2_ttest_annotation$col4 <- 1
turi2_ttest_annotation_out <- select(turi2_ttest_annotation, chemical_iTOL, col1, col2, col3, col4, ttest_color)
turi2_ttest_annotation_out <- turi2_ttest_annotation_out[, c(1, 2, 3, 6, 4, 5)]


turi2_ssmd_annotation$col1 <- 2
turi2_ssmd_annotation$col2 <- 2
turi2_ssmd_annotation$col3 <- 1
turi2_ssmd_annotation$col4 <- 1
turi2_ssmd_annotation_out <- select(turi2_ssmd_annotation, chemical_iTOL, col1, col2, col3, col4, ssmd_color)
turi2_ssmd_annotation_out <- turi2_ssmd_annotation_out[, c(1, 2, 3, 6, 4, 5)]

#write these files out, just need to copy and paste them into a text file for iTOl annotation. the templates can be found on the itol website
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics")

write.table(turi2_ttest_annotation_out, "final_results/turi2_ttest_annotation_out.csv", sep = ",", col.names = F, quote = F, row.names = F)
write.table(turi2_ssmd_annotation_out, "final_results/turi2_ssmd_annotation_out.csv", sep = ",", col.names = F, quote = F, row.names = F)

#extra tabulation and results scripts
#which are categorized differently between the two?
turi2_diffs <- left_join(turi2_ssmd, turi2_ttest, by = "UID")
turi2_diffs <- select(turi2_diffs, UID, ttest_category, ssmd_category)
turi2_diffs <- left_join(turi2_diffs, platemap_id_key_output, by ="UID")
turi2_diffs <- subset(turi2_diffs, (ttest_category != ssmd_category))
#89 compounds
write.csv(turi2_diffs, "final_results/turi2_diffs.csv")

#tabulate for each solvent, define function
#for the ttest
hist_table_ttest <- function(df_result, solv){
  #basically just extract the histogram of the category column, and also output the tabulation
  df_result %>%
  filter(solvent == solv) %>%
  pull(ttest_category) %>%
  hist(.)

  df_result %>%
  filter(solvent == solv) %>%
  select(ttest_category)%>%
  table()
}

hist_table_ttest(turi2_all_ttest_both_category, "DMSO")
hist_table_ttest(turi2_all_ttest_both_category, "MeOH")
hist_table_ttest(turi2_all_ttest_both_category, "H2O")

hist_table_ssmd <- function(df_result, solv){
  #basically just extract the histogram of the category column, and also output the tabulation
  df_result %>%
    filter(solvent == solv) %>%
    pull(ssmd_category) %>%
    hist(.)
  
  df_result %>%
    filter(solvent == solv) %>%
    select(ssmd_category)%>%
    table()
}

hist_table_ssmd(turi2_all_ssmd_both_category, "DMSO")
hist_table_ssmd(turi2_all_ssmd_both_category, "MeOH")
hist_table_ssmd(turi2_all_ssmd_both_category, "H2O")
