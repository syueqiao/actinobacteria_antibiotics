#wroking on making functions
####TO DO#####
#1. make loops for assigning the platemaps and stuff
#2. make neat so that each function can be used for each case -> idk think about this i guess
#3. note that the ggplot uses something different now
#4. generate the dataframes for each strain/solvent


#load packages that are needed
library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(lspline)

#set wd
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/strains_data")
#set this function
`%notin%` <- Negate(`%in%`)

#import platemap
#see supporting files for this document
platemap <- read.csv("../plate_maps.csv", header = TRUE)

#make platemaps for all 3 solvents ###WIP###
solvs <- c("DMSO", "H2O", "MeOH")
for (type in solvs) {
name <- paste("platemap_", type, sep ="")
assign(name, select(platemap, Well, type))
}

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


#################################TO DO
#make inputs for arth######
arth_df_DMSO <- filter(df_metadata, bug == "ArthBac", solvent == "DMSO")
arth_df_H2O <- filter(df_metadata, bug == "ArthBac", solvent == "H2O")

#make inputs for absc######
absc_df <- filter(df_metadata, bug == "absc")
absc_df_DMSO <- filter(absc_df_1 , solvent == "DMSO")
absc_df_H2O <- filter(absc_df_1 , solvent == "H2O")
absc_df_MeOH <- filter(absc_df_1 , solvent == "MeOH")

#make inputs for coryne####
coryn_df <- filter(df_metadata, bug == "coryn")
coryn_df_DMSO <- filter(coryn_df , solvent == "DMSO")
coryn_df_H2O <- filter(coryn_df , solvent == "H2O")
coryn_df_MeOH <- filter(coryn_df , solvent == "MeOH")



###FILTER OUT TIME POINT 1####
#input = whatever dataframe you want to remove the first time point from
remove_first_time <- function(input){
  input_1 <- filter(input, time != "1")
  input_1$time <- as.numeric(input_1$time)
  input_1$time <- input_1$time-1 
  return(input_1)
}


#make dataframe with only the library wells####
### x = whatever final dataframe you want to run the analysis on
analysis_test_DMSO <- function(x){

arth_df_library_DMSO <- filter(x, well %notin% DMSO_filter_out)
#make dataframe with only the control wells
arth_df_bac_DMSO <- filter(x, well %in% DMSO_bac_control_list)
#make dataframe with only the antibiotic wells
arth_df_pos_DMSO <- filter(x, well %in% DMSO_positive_control_list)

#format in usable shape for library
arth_df_library_DMSO <- spread(arth_df_library_DMSO, key = time, value = OD)
#find diff from t8-t1
arth_df_library_DMSO$diff <- arth_df_library_DMSO$"8" - arth_df_library_DMSO$"1"


#find spline from first two time points
#select the values only
arth_df_library_DMSO_splines <- select(arth_df_library_DMSO, c(6:13))

#make time vector of right length
arth_library_time_DMSO <- c(0:7)

#find splines and extract coefficient
knots <- c(3)
arth_lm_DMSO <- apply(arth_df_library_DMSO_splines, 1, function(x) lm(x ~ lspline(arth_library_time_DMSO, knots = knots)))
arth_df_library_DMSO_splines_coeff <- as.data.frame(unlist(lapply(arth_lm_DMSO, function(x) coef(x)[2])))
colnames(arth_df_library_DMSO_splines_coeff) <- c("spline")

#combine with dig dataframe
arth_df_library_DMSO$spline <- arth_df_library_DMSO_splines_coeff

#spread the data for t tests
arth_df_library_DMSO_spread_spline <- spread(arth_df_library_DMSO, key = well, value = spline)
#grab correct number of columns from the table
cols1 <- ncol(arth_df_library_DMSO)-1
cols2 <- ncol(arth_df_library_DMSO_spread_spline)
#adjust to usable format
arth_df_library_DMSO_spread_spline <- select(arth_df_library_DMSO_spread_spline, c(all_of(cols1):all_of(cols2)))
arth_df_library_DMSO_spread_spline <- data.frame(lapply(arth_df_library_DMSO_spread_spline, na.omit))

arth_df_library_DMSO_spread_diff <- spread(arth_df_library_DMSO, key = well, value = diff)
arth_df_library_DMSO_spread_diff <- select(arth_df_library_DMSO_spread_diff, c(all_of(cols1):all_of(cols2)))
arth_df_library_DMSO_spread_diff <- data.frame(lapply(arth_df_library_DMSO_spread_diff, na.omit))

#get vectors for the bac only dataset
#vector for OD
arth_df_bac_DMSO <- spread(arth_df_bac_DMSO, key = time, value = OD)
#filter values above OD 4.5 - example that could work
#arth_df_bac_DMSO <- arth_df_bac_DMSO[arth_df_bac_DMSO$"8" < 4.5]
###################TEST FILTER###########################
arth_df_bac_DMSO$diff <- arth_df_bac_DMSO$"8" - arth_df_bac_DMSO$"1"
arth_df_bac_DMSO_od_vector <- arth_df_bac_DMSO$diff

#vector for spline
arth_df_bac_DMSO_splines <- select(arth_df_bac_DMSO, c(6:13))
arth_bac_time_DMSO <- c(0:7)
arth_bac_lm_DMSO <- apply(arth_df_bac_DMSO_splines, 1, function(x) lm(x ~ lspline(arth_bac_time_DMSO, knots = knots)))
arth_df_bac_DMSO_splines_coeff <- as.data.frame(unlist(lapply(arth_bac_lm_DMSO, function(x) coef(x)[2])))
#combine
arth_df_bac_DMSO$spline <- arth_df_bac_DMSO_splines_coeff
arth_df_bac_DMSO_spline_vector <- arth_df_bac_DMSO$spline

#apply t test over all the columns against the bac-only vectors
#this is for the od diff
arth_ttests_od_DMSO <- lapply(arth_df_library_DMSO_spread_diff, function(x) wilcox.test(x, arth_df_bac_DMSO_od_vector, exact = T)$p.value)
arth_ttests_od_DMSO <- as.data.frame(unlist(arth_ttests_od_DMSO))
colnames(arth_ttests_od_DMSO) <- c("OD diff significance")
#wilcox.test(arth_df_bac_od_vector, arth_df_time_spread_diff$C15, exact = T)

#this is for the splines
#wilcox.test(arth_df_bac_splines_vector, arth_df_time_spread_spline$C15, exact = T)
arth_ttests_spline_DMSO <- lapply(arth_df_library_DMSO_spread_spline, function(x) wilcox.test(x, arth_df_bac_DMSO_spline_vector, exact = T)$p.value)
arth_ttests_spline_DMSO <- as.data.frame(unlist(arth_ttests_spline_DMSO))
colnames(arth_ttests_spline_DMSO) <- c("spline diff significance")

#combine the two together
arth_ttests_combined_DMSO <- cbind(arth_ttests_od_DMSO, arth_ttests_spline_DMSO)

#rownames to column
arth_ttests_combined_DMSO <- rownames_to_column(arth_ttests_combined_DMSO, "well")

#compare?
#make it so that you can compare values
arth_ttests_combined_plus_data_DMSO <- merge(arth_ttests_combined_DMSO, arth_df_library_DMSO, by = "well")
arth_ttests_combined_plus_data_DMSO$od_means <- ave(arth_ttests_combined_plus_data_DMSO$diff, arth_ttests_combined_plus_data_DMSO$well)
arth_ttests_combined_plus_data_DMSO$spline_means <- ave(arth_ttests_combined_plus_data_DMSO$spline, arth_ttests_combined_plus_data_DMSO$well)

#means of control
arth_df_bac_DMSO_od_vector_mean <- mean(arth_df_bac_DMSO_od_vector)
arth_df_bac_DMSO_spline_vector_mean <- mean(arth_df_bac_DMSO_spline_vector)

#compare od diff
arth_ttests_combined_plus_data_DMSO$od_diff <- ifelse(arth_ttests_combined_plus_data_DMSO$od_means > arth_df_bac_DMSO_od_vector_mean, "higher_od", "lower_od")
#compare spline diff
arth_ttests_combined_plus_data_DMSO$spline_diff <- ifelse(arth_ttests_combined_plus_data_DMSO$spline_means > arth_df_bac_DMSO_spline_vector_mean, "higher_spline", "lower_spline")

#determine if differences are statistically significant
arth_ttests_combined_plus_data_DMSO$od_sig <- ifelse(arth_ttests_combined_plus_data_DMSO$"OD diff significance" < 0.05, "sig_od", "ns_od")
arth_ttests_combined_plus_data_DMSO$spline_sig <- ifelse(arth_ttests_combined_plus_data_DMSO$"spline diff significance" < 0.05, "sig_spline", "ns_spline")

#finally, categorize
#this is for things that have no significant difference in either condition from the control 
arth_ttests_combined_plus_data_DMSO$category <- ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "ns_od" & arth_ttests_combined_plus_data_DMSO$spline_sig == "ns_spline", 5, 
                                                       ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "sig_od" & arth_ttests_combined_plus_data_DMSO$od_diff == "lower_od"  & arth_ttests_combined_plus_data_DMSO$spline_sig == "ns_spline", 8, 
                                                              ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "sig_od" & arth_ttests_combined_plus_data_DMSO$od_diff == "higher_od"  & arth_ttests_combined_plus_data_DMSO$spline_sig == "ns_spline", 2, 
                                                                     ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "ns_od" & arth_ttests_combined_plus_data_DMSO$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_DMSO$spline_diff == "higher_spline", 4, 
                                                                            ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "ns_od" & arth_ttests_combined_plus_data_DMSO$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_DMSO$spline_diff == "lower_spline", 6, 
                                                                                   ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "sig_od" & arth_ttests_combined_plus_data_DMSO$od_diff == "higher_od" & arth_ttests_combined_plus_data_DMSO$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_DMSO$spline_diff == "higher_spline", 1,
                                                                                          ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "sig_od" & arth_ttests_combined_plus_data_DMSO$od_diff == "lower_od" & arth_ttests_combined_plus_data_DMSO$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_DMSO$spline_diff == "higher_spline", 7, 
                                                                                                 ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "sig_od" & arth_ttests_combined_plus_data_DMSO$od_diff == "higher_od" & arth_ttests_combined_plus_data_DMSO$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_DMSO$spline_diff == "lower_spline", 3,
                                                                                                        ifelse(arth_ttests_combined_plus_data_DMSO$od_sig == "sig_od" & arth_ttests_combined_plus_data_DMSO$od_diff == "lower_od" & arth_ttests_combined_plus_data_DMSO$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_DMSO$spline_diff == "lower_spline", 9, 0)))))))))



#quick qc
hist(arth_ttests_combined_plus_data_DMSO$category, breaks = 12)
table(arth_ttests_combined_plus_data_DMSO$category)
arth_ttests_histo <- c(arth_ttests_combined_plus_data_DMSO$category)

#make nicer on the eyes
ttest_output_arth_DMSO <- select(arth_ttests_combined_plus_data_DMSO, well, category)
ttest_output_arth_DMSO <- ttest_output_arth_DMSO[!duplicated(ttest_output_arth_DMSO), ]
arth_DMSO_library <- platemap_DMSO
names(arth_DMSO_library)[names(arth_DMSO_library) == 'Well'] <- 'well'

#merge output with library files
ttest_output_arth_DMSO <<- merge(arth_DMSO_library, ttest_output_arth_DMSO, by = "well" )
return(ttest_output_arth_DMSO)
}
analysis_test_H2O <- function(x){
  
  arth_df_library_H2O <- filter(x, well %notin% H2O_filter_out)
  #make dataframe with only the control wells
  arth_df_bac_H2O <- filter(x, well %in% H2O_bac_control_list)
  #make dataframe with only the antibiotic wells
  arth_df_pos_H2O <- filter(x, well %in% H2O_positive_control_list)
  
  #format in usable shape for library
  arth_df_library_H2O <- spread(arth_df_library_H2O, key = time, value = OD)
  #find diff from t8-t1
  arth_df_library_H2O$diff <- arth_df_library_H2O$"8" - arth_df_library_H2O$"1"
  
  
  #find spline from first two time points
  #select the values only
  arth_df_library_H2O_splines <- select(arth_df_library_H2O, c(6:13))
  
  #make time vector of right length
  arth_library_time_H2O <- c(0:7)
  
  #find splines and extract coefficient
  knots <- c(3)
  arth_lm_H2O <- apply(arth_df_library_H2O_splines, 1, function(x) lm(x ~ lspline(arth_library_time_H2O, knots = knots)))
  arth_df_library_H2O_splines_coeff <- as.data.frame(unlist(lapply(arth_lm_H2O, function(x) coef(x)[2])))
  colnames(arth_df_library_H2O_splines_coeff) <- c("spline")
  
  #combine with dig dataframe
  arth_df_library_H2O$spline <- arth_df_library_H2O_splines_coeff
  
  #spread the data for t tests
  arth_df_library_H2O_spread_spline <- spread(arth_df_library_H2O, key = well, value = spline)
  #grab correct number of columns from the table
  cols1 <- ncol(arth_df_library_H2O)-1
  cols2 <- ncol(arth_df_library_H2O_spread_spline)
  #adjust to usable format
  arth_df_library_H2O_spread_spline <- select(arth_df_library_H2O_spread_spline, c(all_of(cols1):all_of(cols2)))
  arth_df_library_H2O_spread_spline <- data.frame(lapply(arth_df_library_H2O_spread_spline, na.omit))
  
  arth_df_library_H2O_spread_diff <- spread(arth_df_library_H2O, key = well, value = diff)
  arth_df_library_H2O_spread_diff <- select(arth_df_library_H2O_spread_diff, c(all_of(cols1):all_of(cols2)))
  arth_df_library_H2O_spread_diff <- data.frame(lapply(arth_df_library_H2O_spread_diff, na.omit))
  
  #get vectors for the bac only dataset
  #vector for OD
  arth_df_bac_H2O <- spread(arth_df_bac_H2O, key = time, value = OD)
  #filter values above OD 4.5 - example that could work
  #arth_df_bac_H2O <- arth_df_bac_H2O[arth_df_bac_H2O$"8" < 4.5]
  ###################TEST FILTER###########################
  arth_df_bac_H2O$diff <- arth_df_bac_H2O$"8" - arth_df_bac_H2O$"1"
  arth_df_bac_H2O_od_vector <- arth_df_bac_H2O$diff
  
  #vector for spline
  arth_df_bac_H2O_splines <- select(arth_df_bac_H2O, c(6:13))
  arth_bac_time_H2O <- c(0:7)
  arth_bac_lm_H2O <- apply(arth_df_bac_H2O_splines, 1, function(x) lm(x ~ lspline(arth_bac_time_H2O, knots = knots)))
  arth_df_bac_H2O_splines_coeff <- as.data.frame(unlist(lapply(arth_bac_lm_H2O, function(x) coef(x)[2])))
  #combine
  arth_df_bac_H2O$spline <- arth_df_bac_H2O_splines_coeff
  arth_df_bac_H2O_spline_vector <- arth_df_bac_H2O$spline
  
  #apply t test over all the columns against the bac-only vectors
  #this is for the od diff
  arth_ttests_od_H2O <- lapply(arth_df_library_H2O_spread_diff, function(x) wilcox.test(x, arth_df_bac_H2O_od_vector, exact = T)$p.value)
  arth_ttests_od_H2O <- as.data.frame(unlist(arth_ttests_od_H2O))
  colnames(arth_ttests_od_H2O) <- c("OD diff significance")
  #wilcox.test(arth_df_bac_od_vector, arth_df_time_spread_diff$C15, exact = T)
  
  #this is for the splines
  #wilcox.test(arth_df_bac_splines_vector, arth_df_time_spread_spline$C15, exact = T)
  arth_ttests_spline_H2O <- lapply(arth_df_library_H2O_spread_spline, function(x) wilcox.test(x, arth_df_bac_H2O_spline_vector, exact = T)$p.value)
  arth_ttests_spline_H2O <- as.data.frame(unlist(arth_ttests_spline_H2O))
  colnames(arth_ttests_spline_H2O) <- c("spline diff significance")
  
  #combine the two together
  arth_ttests_combined_H2O <- cbind(arth_ttests_od_H2O, arth_ttests_spline_H2O)
  
  #rownames to column
  arth_ttests_combined_H2O <- rownames_to_column(arth_ttests_combined_H2O, "well")
  
  #compare?
  #make it so that you can compare values
  arth_ttests_combined_plus_data_H2O <- merge(arth_ttests_combined_H2O, arth_df_library_H2O, by = "well")
  arth_ttests_combined_plus_data_H2O$od_means <- ave(arth_ttests_combined_plus_data_H2O$diff, arth_ttests_combined_plus_data_H2O$well)
  arth_ttests_combined_plus_data_H2O$spline_means <- ave(arth_ttests_combined_plus_data_H2O$spline, arth_ttests_combined_plus_data_H2O$well)
  
  #means of control
  arth_df_bac_H2O_od_vector_mean <- mean(arth_df_bac_H2O_od_vector)
  arth_df_bac_H2O_spline_vector_mean <- mean(arth_df_bac_H2O_spline_vector)
  
  #compare od diff
  arth_ttests_combined_plus_data_H2O$od_diff <- ifelse(arth_ttests_combined_plus_data_H2O$od_means > arth_df_bac_H2O_od_vector_mean, "higher_od", "lower_od")
  #compare spline diff
  arth_ttests_combined_plus_data_H2O$spline_diff <- ifelse(arth_ttests_combined_plus_data_H2O$spline_means > arth_df_bac_H2O_spline_vector_mean, "higher_spline", "lower_spline")
  
  #determine if differences are statistically significant
  arth_ttests_combined_plus_data_H2O$od_sig <- ifelse(arth_ttests_combined_plus_data_H2O$"OD diff significance" < 0.05, "sig_od", "ns_od")
  arth_ttests_combined_plus_data_H2O$spline_sig <- ifelse(arth_ttests_combined_plus_data_H2O$"spline diff significance" < 0.05, "sig_spline", "ns_spline")
  
  #finally, categorize
  #this is for things that have no significant difference in either condition from the control 
  arth_ttests_combined_plus_data_H2O$category <- ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "ns_od" & arth_ttests_combined_plus_data_H2O$spline_sig == "ns_spline", 5, 
                                                         ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "sig_od" & arth_ttests_combined_plus_data_H2O$od_diff == "lower_od"  & arth_ttests_combined_plus_data_H2O$spline_sig == "ns_spline", 8, 
                                                                ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "sig_od" & arth_ttests_combined_plus_data_H2O$od_diff == "higher_od"  & arth_ttests_combined_plus_data_H2O$spline_sig == "ns_spline", 2, 
                                                                       ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "ns_od" & arth_ttests_combined_plus_data_H2O$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_H2O$spline_diff == "higher_spline", 4, 
                                                                              ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "ns_od" & arth_ttests_combined_plus_data_H2O$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_H2O$spline_diff == "lower_spline", 6, 
                                                                                     ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "sig_od" & arth_ttests_combined_plus_data_H2O$od_diff == "higher_od" & arth_ttests_combined_plus_data_H2O$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_H2O$spline_diff == "higher_spline", 1,
                                                                                            ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "sig_od" & arth_ttests_combined_plus_data_H2O$od_diff == "lower_od" & arth_ttests_combined_plus_data_H2O$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_H2O$spline_diff == "higher_spline", 7, 
                                                                                                   ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "sig_od" & arth_ttests_combined_plus_data_H2O$od_diff == "higher_od" & arth_ttests_combined_plus_data_H2O$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_H2O$spline_diff == "lower_spline", 3,
                                                                                                          ifelse(arth_ttests_combined_plus_data_H2O$od_sig == "sig_od" & arth_ttests_combined_plus_data_H2O$od_diff == "lower_od" & arth_ttests_combined_plus_data_H2O$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_H2O$spline_diff == "lower_spline", 9, 0)))))))))
  
  
  
  #quick qc
  hist(arth_ttests_combined_plus_data_H2O$category, breaks = 12)
  table(arth_ttests_combined_plus_data_H2O$category)
  arth_ttests_histo <- c(arth_ttests_combined_plus_data_H2O$category)
  
  #make nicer on the eyes
  ttest_output_arth_H2O <- select(arth_ttests_combined_plus_data_H2O, well, category)
  ttest_output_arth_H2O <- ttest_output_arth_H2O[!duplicated(ttest_output_arth_H2O), ]
  arth_H2O_library <- platemap_H2O
  names(arth_H2O_library)[names(arth_H2O_library) == 'Well'] <- 'well'
  
  #merge output with library files
  ttest_output_arth_H2O <<- merge(arth_H2O_library, ttest_output_arth_H2O, by = "well" )
  return(ttest_output_arth_H2O)
}
analysis_test_MeOH <- function(x){
  
  arth_df_library_MeOH <- filter(x, well %notin% MeOH_filter_out)
  #make dataframe with only the control wells
  arth_df_bac_MeOH <- filter(x, well %in% MeOH_bac_control_list)
  #make dataframe with only the antibiotic wells
  arth_df_pos_MeOH <- filter(x, well %in% MeOH_positive_control_list)
  
  #format in usable shape for library
  arth_df_library_MeOH <- spread(arth_df_library_MeOH, key = time, value = OD)
  #find diff from t8-t1
  arth_df_library_MeOH$diff <- arth_df_library_MeOH$"8" - arth_df_library_MeOH$"1"
  
  
  #find spline from first two time points
  #select the values only
  arth_df_library_MeOH_splines <- select(arth_df_library_MeOH, c(6:13))
  
  #make time vector of right length
  arth_library_time_MeOH <- c(0:7)
  
  #find splines and extract coefficient
  knots <- c(3)
  arth_lm_MeOH <- apply(arth_df_library_MeOH_splines, 1, function(x) lm(x ~ lspline(arth_library_time_MeOH, knots = knots)))
  arth_df_library_MeOH_splines_coeff <- as.data.frame(unlist(lapply(arth_lm_MeOH, function(x) coef(x)[2])))
  colnames(arth_df_library_MeOH_splines_coeff) <- c("spline")
  
  #combine with dig dataframe
  arth_df_library_MeOH$spline <- arth_df_library_MeOH_splines_coeff
  
  #spread the data for t tests
  arth_df_library_MeOH_spread_spline <- spread(arth_df_library_MeOH, key = well, value = spline)
  #grab correct number of columns from the table
  cols1 <- ncol(arth_df_library_MeOH)-1
  cols2 <- ncol(arth_df_library_MeOH_spread_spline)
  #adjust to usable format
  arth_df_library_MeOH_spread_spline <- select(arth_df_library_MeOH_spread_spline, c(all_of(cols1):all_of(cols2)))
  arth_df_library_MeOH_spread_spline <- data.frame(lapply(arth_df_library_MeOH_spread_spline, na.omit))
  
  arth_df_library_MeOH_spread_diff <- spread(arth_df_library_MeOH, key = well, value = diff)
  arth_df_library_MeOH_spread_diff <- select(arth_df_library_MeOH_spread_diff, c(all_of(cols1):all_of(cols2)))
  arth_df_library_MeOH_spread_diff <- data.frame(lapply(arth_df_library_MeOH_spread_diff, na.omit))
  
  #get vectors for the bac only dataset
  #vector for OD
  arth_df_bac_MeOH <- spread(arth_df_bac_MeOH, key = time, value = OD)
  #filter values above OD 4.5 - example that could work
  #arth_df_bac_MeOH <- arth_df_bac_MeOH[arth_df_bac_MeOH$"8" < 4.5]
  ###################TEST FILTER###########################
  arth_df_bac_MeOH$diff <- arth_df_bac_MeOH$"8" - arth_df_bac_MeOH$"1"
  arth_df_bac_MeOH_od_vector <- arth_df_bac_MeOH$diff
  
  #vector for spline
  arth_df_bac_MeOH_splines <- select(arth_df_bac_MeOH, c(6:13))
  arth_bac_time_MeOH <- c(0:7)
  arth_bac_lm_MeOH <- apply(arth_df_bac_MeOH_splines, 1, function(x) lm(x ~ lspline(arth_bac_time_MeOH, knots = knots)))
  arth_df_bac_MeOH_splines_coeff <- as.data.frame(unlist(lapply(arth_bac_lm_MeOH, function(x) coef(x)[2])))
  #combine
  arth_df_bac_MeOH$spline <- arth_df_bac_MeOH_splines_coeff
  arth_df_bac_MeOH_spline_vector <- arth_df_bac_MeOH$spline
  
  #apply t test over all the columns against the bac-only vectors
  #this is for the od diff
  arth_ttests_od_MeOH <- lapply(arth_df_library_MeOH_spread_diff, function(x) wilcox.test(x, arth_df_bac_MeOH_od_vector, exact = T)$p.value)
  arth_ttests_od_MeOH <- as.data.frame(unlist(arth_ttests_od_MeOH))
  colnames(arth_ttests_od_MeOH) <- c("OD diff significance")
  #wilcox.test(arth_df_bac_od_vector, arth_df_time_spread_diff$C15, exact = T)
  
  #this is for the splines
  #wilcox.test(arth_df_bac_splines_vector, arth_df_time_spread_spline$C15, exact = T)
  arth_ttests_spline_MeOH <- lapply(arth_df_library_MeOH_spread_spline, function(x) wilcox.test(x, arth_df_bac_MeOH_spline_vector, exact = T)$p.value)
  arth_ttests_spline_MeOH <- as.data.frame(unlist(arth_ttests_spline_MeOH))
  colnames(arth_ttests_spline_MeOH) <- c("spline diff significance")
  
  #combine the two together
  arth_ttests_combined_MeOH <- cbind(arth_ttests_od_MeOH, arth_ttests_spline_MeOH)
  
  #rownames to column
  arth_ttests_combined_MeOH <- rownames_to_column(arth_ttests_combined_MeOH, "well")
  
  #compare?
  #make it so that you can compare values
  arth_ttests_combined_plus_data_MeOH <- merge(arth_ttests_combined_MeOH, arth_df_library_MeOH, by = "well")
  arth_ttests_combined_plus_data_MeOH$od_means <- ave(arth_ttests_combined_plus_data_MeOH$diff, arth_ttests_combined_plus_data_MeOH$well)
  arth_ttests_combined_plus_data_MeOH$spline_means <- ave(arth_ttests_combined_plus_data_MeOH$spline, arth_ttests_combined_plus_data_MeOH$well)
  
  #means of control
  arth_df_bac_MeOH_od_vector_mean <- mean(arth_df_bac_MeOH_od_vector)
  arth_df_bac_MeOH_spline_vector_mean <- mean(arth_df_bac_MeOH_spline_vector)
  
  #compare od diff
  arth_ttests_combined_plus_data_MeOH$od_diff <- ifelse(arth_ttests_combined_plus_data_MeOH$od_means > arth_df_bac_MeOH_od_vector_mean, "higher_od", "lower_od")
  #compare spline diff
  arth_ttests_combined_plus_data_MeOH$spline_diff <- ifelse(arth_ttests_combined_plus_data_MeOH$spline_means > arth_df_bac_MeOH_spline_vector_mean, "higher_spline", "lower_spline")
  
  #determine if differences are statistically significant
  arth_ttests_combined_plus_data_MeOH$od_sig <- ifelse(arth_ttests_combined_plus_data_MeOH$"OD diff significance" < 0.05, "sig_od", "ns_od")
  arth_ttests_combined_plus_data_MeOH$spline_sig <- ifelse(arth_ttests_combined_plus_data_MeOH$"spline diff significance" < 0.05, "sig_spline", "ns_spline")
  
  #finally, categorize
  #this is for things that have no significant difference in either condition from the control 
  arth_ttests_combined_plus_data_MeOH$category <- ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "ns_od" & arth_ttests_combined_plus_data_MeOH$spline_sig == "ns_spline", 5, 
                                                        ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "sig_od" & arth_ttests_combined_plus_data_MeOH$od_diff == "lower_od"  & arth_ttests_combined_plus_data_MeOH$spline_sig == "ns_spline", 8, 
                                                               ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "sig_od" & arth_ttests_combined_plus_data_MeOH$od_diff == "higher_od"  & arth_ttests_combined_plus_data_MeOH$spline_sig == "ns_spline", 2, 
                                                                      ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "ns_od" & arth_ttests_combined_plus_data_MeOH$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_MeOH$spline_diff == "higher_spline", 4, 
                                                                             ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "ns_od" & arth_ttests_combined_plus_data_MeOH$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_MeOH$spline_diff == "lower_spline", 6, 
                                                                                    ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "sig_od" & arth_ttests_combined_plus_data_MeOH$od_diff == "higher_od" & arth_ttests_combined_plus_data_MeOH$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_MeOH$spline_diff == "higher_spline", 1,
                                                                                           ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "sig_od" & arth_ttests_combined_plus_data_MeOH$od_diff == "lower_od" & arth_ttests_combined_plus_data_MeOH$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_MeOH$spline_diff == "higher_spline", 7, 
                                                                                                  ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "sig_od" & arth_ttests_combined_plus_data_MeOH$od_diff == "higher_od" & arth_ttests_combined_plus_data_MeOH$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_MeOH$spline_diff == "lower_spline", 3,
                                                                                                         ifelse(arth_ttests_combined_plus_data_MeOH$od_sig == "sig_od" & arth_ttests_combined_plus_data_MeOH$od_diff == "lower_od" & arth_ttests_combined_plus_data_MeOH$spline_sig == "sig_spline" & arth_ttests_combined_plus_data_MeOH$spline_diff == "lower_spline", 9, 0)))))))))
  
  
  
  #quick qc
  hist(arth_ttests_combined_plus_data_MeOH$category, breaks = 12)
  table(arth_ttests_combined_plus_data_MeOH$category)
  arth_ttests_histo <- c(arth_ttests_combined_plus_data_MeOH$category)
  
  #make nicer on the eyes
  ttest_output_arth_MeOH <- select(arth_ttests_combined_plus_data_MeOH, well, category)
  ttest_output_arth_MeOH <- ttest_output_arth_MeOH[!duplicated(ttest_output_arth_MeOH), ]
  arth_MeOH_library <- platemap_MeOH
  names(arth_MeOH_library)[names(arth_MeOH_library) == 'Well'] <- 'well'
  
  #merge output with library files
  ttest_output_arth_MeOH <<- merge(arth_MeOH_library, ttest_output_arth_MeOH, by = "well" )
  return(ttest_output_arth_MeOH)
}
#don't forget to assign output to whatever you run through the function

#can check the numerical values using 
cat_numbers <- function(cat){
  return(table(cat$category))
}
#export data
path <- 'C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/category_outputs'
write.csv(ttest_output_arth_DMSO, file.path(path, paste("arth_DMSO_categories.csv", sep = ""), row.names = F))


#plot prelim figure#######

#input = output of the earlier function, don't forget to change the name of the graph title when exporting
ggplot(absc_DMSO, aes(x=category, fill=..x..)) +
  geom_bar() +
  xlab("Category") +
  ylab("Count") +
  ggtitle("Effect on Arthrobacter growth of DMSO-solvent compounds") +
  scale_x_continuous(breaks = 1:9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = '#5AAA46', high = '#825CA6', breaks = c(1:9), labels = c(1:9), name = "Categories")



