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
###############make inputs for arth#################
arth_df_DMSO <- filter(df_metadata, bug == "ArthBac", solvent == "DMSO")
arth_df_H2O <- filter(df_metadata, bug == "ArthBac", solvent == "H2O")
arth_df_MeOH <- filter(df_metadata, bug == "ArthBac", solvent == "MeOH")
#look at the graphs manually, and determine data to be filtered out
#for example, for arthrobacter_DMSO, I want to filter out bacteria-only that are above OD 3, can just remove single replicates for this, since we consider all reps part of one sample
#I also want to filter out library compounds that have an OD above 5
#very roundabout way but couldn't think of a better way to do it
arth_df_library_DMSO <- filter(arth_df_DMSO, well %notin% DMSO_filter_out)
arth_df_library_DMSO_spread <- spread(arth_df_library_DMSO, key = time, value = OD)
arth_df_library_out <- filter(arth_df_library_DMSO_spread, arth_df_library_DMSO_spread$"8" > 5)
arth_df_library_out <- arth_df_library_out$well
arth_df_library_DMSO <- filter(arth_df_library_DMSO_spread, well %notin% arth_df_library_out)

#make dataframe with only the control wells
arth_df_bac_DMSO <- filter(arth_df_DMSO, well %in% DMSO_bac_control_list)
arth_df_bac_DMSO_spread <- spread(arth_df_bac_DMSO, key = time, value = OD)
arth_df_bac_DMSO <- filter(arth_df_bac_DMSO_spread, arth_df_bac_DMSO_spread$"8" < 3)

#make dataframe with only the antibiot ic wells
arth_df_pos_DMSO <- filter(arth_df_DMSO, well %in% DMSO_positive_control_list)

#make similar dataframe but for MeOH
arth_df_library_MeOH <- filter(arth_df_MeOH, well %notin% MeOH_filter_out)
arth_df_library_MeOH_spread <- spread(arth_df_library_MeOH, key = time, value = OD)
arth_df_library_out <- filter(arth_df_library_MeOH_spread, arth_df_library_MeOH_spread$"8" > 5)
arth_df_library_out <- arth_df_library_out$well
arth_df_library_MeOH <- filter(arth_df_library_MeOH_spread, well %notin% arth_df_library_out)

#make dataframe with only the control wells
arth_df_bac_MeOH <- filter(arth_df_MeOH, well %in% MeOH_bac_control_list)
arth_df_bac_MeOH_spread <- spread(arth_df_bac_MeOH, key = time, value = OD)
arth_df_bac_MeOH <- filter(arth_df_bac_MeOH_spread, arth_df_bac_MeOH_spread$"8" < 2 & arth_df_bac_MeOH_spread$"8" > 1)

#make dataframe with only the antibiot ic wells
arth_df_pos_MeOH <- filter(arth_df_MeOH, well %in% MeOH_positive_control_list)

#make for H2O
arth_df_library_H2O <- filter(arth_dsf_H2O, well %notin% H2O_filter_out)
arth_df_library_H2O_spread <- spread(arth_df_library_H2O, key = time, value = OD)
arth_df_library_out <- filter(arth_df_library_H2O_spread, arth_df_library_H2O_spread$"8" > 5)
arth_df_library_out <- arth_df_library_out$well
arth_df_library_H2O <- filter(arth_df_library_H2O_spread, well %notin% arth_df_library_out)

#make dataframe with only the control wells
arth_df_bac_H2O <- filter(arth_df_H2O, well %in% H2O_bac_control_list)
arth_df_bac_H2O_spread <- spread(arth_df_bac_H2O, key = time, value = OD)
arth_df_bac_H2O <- filter(arth_df_bac_H2O_spread, arth_df_bac_H2O_spread$"8" < 3 & arth_df_bac_H2O_spread$"8" > 1)

#make dataframe with only the antibiot ic wells
arth_df_pos_H2O <- filter(arth_df_H2O, well %in% H2O_positive_control_list)

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

coryn_df_library_DMSO <- filter(coryn_df_DMSO, well %notin% DMSO_filter_out)
coryn_df_library_DMSO_spread <- spread(coryn_df_library_DMSO, key = time, value = OD)
coryn_df_library_out <- filter(coryn_df_library_DMSO_spread, coryn_df_library_DMSO_spread$"8" > 3.5 | coryn_df_library_DMSO_spread$"1" > 0.5)
coryn_df_library_out <- coryn_df_library_out$well
coryn_df_library_DMSO <- filter(coryn_df_library_DMSO_spread, well %notin% coryn_df_library_out)

#make dataframe with only the control wells
coryn_df_bac_DMSO <- filter(coryn_df_DMSO, well %in% DMSO_bac_control_list)
coryn_df_bac_DMSO_spread <- spread(coryn_df_bac_DMSO, key = time, value = OD)
coryn_df_bac_DMSO <- filter(coryn_df_bac_DMSO_spread, coryn_df_bac_DMSO_spread$"8" < 1.5)

#make dataframe with only the antibiot ic wells
coryn_df_pos_DMSO <- filter(coryn_df_DMSO, well %in% DMSO_positive_control_list)

#make similar dataframe but for MeOH
coryn_df_library_MeOH <- filter(coryn_df_MeOH, well %notin% MeOH_filter_out)
coryn_df_library_MeOH_spread <- spread(coryn_df_library_MeOH, key = time, value = OD)
coryn_df_library_out <- filter(coryn_df_library_MeOH_spread, coryn_df_library_MeOH_spread$"8" > 3.25)
coryn_df_library_out <- coryn_df_library_out$well
coryn_df_library_MeOH <- filter(coryn_df_library_MeOH_spread, well %notin% coryn_df_library_out)

#make dataframe with only the control wells
coryn_df_bac_MeOH <- filter(coryn_df_MeOH, well %in% MeOH_bac_control_list)
coryn_df_bac_MeOH_spread <- spread(coryn_df_bac_MeOH, key = time, value = OD)
coryn_df_bac_MeOH <- filter(coryn_df_bac_MeOH_spread, coryn_df_bac_MeOH_spread$"8" < 2.5 & coryn_df_bac_MeOH_spread$"8" > 1.4)

#make dataframe with only the antibiot ic wells
coryn_df_pos_MeOH <- filter(coryn_df_MeOH, well %in% MeOH_positive_control_list)

#make for H2O
coryn_df_library_H2O <- filter(coryn_df_H2O, well %notin% H2O_filter_out)
coryn_df_library_H2O_spread <- spread(coryn_df_library_H2O, key = time, value = OD)
coryn_df_library_out <- filter(coryn_df_library_H2O_spread, coryn_df_library_H2O_spread$"8" > 3)
coryn_df_library_out <- coryn_df_library_out$well
coryn_df_library_H2O <- filter(coryn_df_library_H2O_spread, well %notin% coryn_df_library_out)

#make dataframe with only the control wells
coryn_df_bac_H2O <- filter(coryn_df_H2O, well %in% H2O_bac_control_list)
coryn_df_bac_H2O_spread <- spread(coryn_df_bac_H2O, key = time, value = OD)
coryn_df_bac_H2O <- filter(coryn_df_bac_H2O_spread, coryn_df_bac_H2O_spread$"8" < 2 & coryn_df_bac_H2O_spread$"8" > 1)

#make dataframe with only the antibiot ic wells
coryn_df_pos_H2O <- filter(coryn_df_H2O, well %in% H2O_positive_control_list)

#make inputs for rhodo####
rhodo_df <- filter(df_metadata, bug == "rhodo")
rhodo_df_DMSO <- filter(rhodo_df , solvent == "DMSO")
rhodo_df_H2O <- filter(rhodo_df , solvent == "H2O")
rhodo_df_MeOH <- filter(rhodo_df , solvent == "MeOH")

#make inputs for brevi#####
brevi_df <- filter(df_metadata, bug == "brevi")
brevi_df_DMSO <- filter(brevi_df , solvent == "DMSO")
brevi_df_H2O <- filter(brevi_df , solvent == "H2O")
brevi_df_MeOH <- filter(brevi_df , solvent == "MeOH")

brevi_df_library_DMSO <- filter(brevi_df_DMSO, well %notin% DMSO_filter_out)
brevi_df_library_DMSO <- remove_first_time(brevi_df_library_DMSO)
brevi_df_library_DMSO_spread <- spread(brevi_df_library_DMSO, key = time, value = OD)
#this one is different -> filters OUT wells that meet these conditions
brevi_df_library_out <- filter(brevi_df_library_DMSO_spread, brevi_df_library_DMSO_spread$"8" > 5)
brevi_df_library_out <- brevi_df_library_out$well
brevi_df_library_DMSO <- filter(brevi_df_library_DMSO_spread, well %notin% brevi_df_library_out)

#make dataframe with only the control wells
brevi_df_bac_DMSO <- filter(brevi_df_DMSO, well %in% DMSO_bac_control_list)
brevi_df_bac_DMSO <- remove_first_time(brevi_df_bac_DMSO)
brevi_df_bac_DMSO_spread <- spread(brevi_df_bac_DMSO, key = time, value = OD)
brevi_df_bac_DMSO <- filter(brevi_df_bac_DMSO_spread, brevi_df_bac_DMSO_spread$"8" > 3)

#make dataframe with only the antibiot ic wells
brevi_df_pos_DMSO <- filter(brevi_df_DMSO, well %in% DMSO_positive_control_list)

#make similar dataframe but for MeOH
brevi_df_library_MeOH <- filter(brevi_df_MeOH, well %notin% MeOH_filter_out)
brevi_df_library_MeOH <- remove_first_time(brevi_df_library_MeOH)
brevi_df_library_MeOH_spread <- spread(brevi_df_library_MeOH, key = time, value = OD)
brevi_df_library_out <- filter(brevi_df_library_MeOH_spread, brevi_df_library_MeOH_spread$"8" > 4.5)
brevi_df_library_out <- brevi_df_library_out$well
brevi_df_library_MeOH <- filter(brevi_df_library_MeOH_spread, well %notin% brevi_df_library_out)

#make dataframe with only the control wells
brevi_df_bac_MeOH <- filter(brevi_df_MeOH, well %in% MeOH_bac_control_list)
brevi_df_bac_MeOH <- remove_first_time(brevi_df_bac_MeOH)
brevi_df_bac_MeOH_spread <- spread(brevi_df_bac_MeOH, key = time, value = OD)
brevi_df_bac_MeOH <- filter(brevi_df_bac_MeOH_spread, brevi_df_bac_MeOH_spread$"8" > 0.5 & brevi_df_bac_MeOH_spread$"8" < 4)

#make dataframe with only the antibiot ic wells
brevi_df_pos_MeOH <- filter(brevi_df_MeOH, well %in% MeOH_positive_control_list)

#make for H2O
brevi_df_library_H2O <- filter(brevi_df_H2O, well %notin% H2O_filter_out)
brevi_df_library_H2O <- remove_first_time(brevi_df_library_H2O)
brevi_df_library_H2O_spread <- spread(brevi_df_library_H2O, key = time, value = OD)
brevi_df_library_out <- filter(brevi_df_library_H2O_spread, brevi_df_library_H2O_spread$"8" > 5 | brevi_df_library_H2O_spread$"1" > 1)
brevi_df_library_out <- brevi_df_library_out$well
brevi_df_library_H2O <- filter(brevi_df_library_H2O_spread, well %notin% brevi_df_library_out)

#make dataframe with only the control wells
brevi_df_bac_H2O <- filter(brevi_df_H2O, well %in% H2O_bac_control_list)
brevi_df_bac_H2O <- remove_first_time(brevi_df_bac_H2O)
brevi_df_bac_H2O_spread <- spread(brevi_df_bac_H2O, key = time, value = OD)
brevi_df_bac_H2O <- filter(brevi_df_bac_H2O_spread, brevi_df_bac_H2O_spread$"8" < 4 & brevi_df_bac_H2O_spread$"8" > 2.5)

#make dataframe with only the antibiot ic wells
brevi_df_pos_H2O <- filter(brevi_df_H2O, well %in% H2O_positive_control_list)

#make inputs for smeg3#####
smeg3_df <- filter(df_metadata, bug == "smegmatis3")
smeg3_df_DMSO <- filter(smeg3_df , solvent == "DMSO")
smeg3_df_H2O <- filter(smeg3_df , solvent == "H2O")
smeg3_df_MeOH <- filter(smeg3_df , solvent == "MeOH")

smeg3_df_library_DMSO <- filter(smeg3_df_DMSO, well %notin% DMSO_filter_out)
smeg3_df_library_DMSO <- remove_first_time(smeg3_df_library_DMSO)
smeg3_df_library_DMSO <- remove_first_time(smeg3_df_library_DMSO)
smeg3_df_library_DMSO_spread <- spread(smeg3_df_library_DMSO, key = time, value = OD)
smeg3_df_library_out <- filter(smeg3_df_library_DMSO_spread, smeg3_df_library_DMSO_spread$"8" > 1.5 )
smeg3_df_library_out <- smeg3_df_library_out$well
smeg3_df_library_DMSO <- filter(smeg3_df_library_DMSO_spread, well %notin% smeg3_df_library_out)

#make dataframe with only the control wells
smeg3_df_bac_DMSO <- filter(smeg3_df_DMSO, well %in% DMSO_bac_control_list)
smeg3_df_bac_DMSO <- remove_first_time(smeg3_df_bac_DMSO)
smeg3_df_bac_DMSO <- remove_first_time(smeg3_df_bac_DMSO)
smeg3_df_bac_DMSO_spread <- spread(smeg3_df_bac_DMSO, key = time, value = OD)
smeg3_df_bac_DMSO <- filter(smeg3_df_bac_DMSO_spread, smeg3_df_bac_DMSO_spread$"8" < 0.9)

#make dataframe with only the antibiot ic wells
smeg3_df_pos_DMSO <- filter(smeg3_df_DMSO, well %in% DMSO_positive_control_list)

#make similar dataframe but for MeOH
smeg3_df_library_MeOH <- filter(smeg3_df_MeOH, well %notin% MeOH_filter_out)
smeg3_df_library_MeOH <- remove_first_time(smeg3_df_library_MeOH)
smeg3_df_library_MeOH <- remove_first_time(smeg3_df_library_MeOH)
smeg3_df_library_MeOH_spread <- spread(smeg3_df_library_MeOH, key = time, value = OD)
smeg3_df_library_out <- filter(smeg3_df_library_MeOH_spread, smeg3_df_library_MeOH_spread$"8" > 1.25 |smeg3_df_library_MeOH_spread$"1" > 0.25)
smeg3_df_library_out <- smeg3_df_library_out$well
smeg3_df_library_MeOH <- filter(smeg3_df_library_MeOH_spread, well %notin% smeg3_df_library_out)

#make dataframe with only the control wells
smeg3_df_bac_MeOH <- filter(smeg3_df_MeOH, well %in% MeOH_bac_control_list)
smeg3_df_bac_MeOH <- remove_first_time(smeg3_df_bac_MeOH)
smeg3_df_bac_MeOH <- remove_first_time(smeg3_df_bac_MeOH)
smeg3_df_bac_MeOH_spread <- spread(smeg3_df_bac_MeOH, key = time, value = OD)
smeg3_df_bac_MeOH <- filter(smeg3_df_bac_MeOH_spread, smeg3_df_bac_MeOH_spread$"8" < 1)

#make dataframe with only the antibiot ic wells
smeg3_df_pos_MeOH <- filter(smeg3_df_MeOH, well %in% MeOH_positive_control_list)

#make for H2O
smeg3_df_library_H2O <- filter(smeg3_df_H2O, well %notin% H2O_filter_out)
smeg3_df_library_H2O <- remove_first_time(smeg3_df_library_H2O)
smeg3_df_library_H2O <- remove_first_time(smeg3_df_library_H2O)
smeg3_df_library_H2O_spread <- spread(smeg3_df_library_H2O, key = time, value = OD)
smeg3_df_library_out <- filter(smeg3_df_library_H2O_spread, smeg3_df_library_H2O_spread$"8" > 1.25 | smeg3_df_library_H2O_spread$"1" > 0.3)
smeg3_df_library_out <- smeg3_df_library_out$well
smeg3_df_library_H2O <- filter(smeg3_df_library_H2O_spread, well %notin% smeg3_df_library_out)

#make dataframe with only the control wells
smeg3_df_bac_H2O <- filter(smeg3_df_H2O, well %in% H2O_bac_control_list)
smeg3_df_bac_H2O <- remove_first_time(smeg3_df_bac_H2O)
smeg3_df_bac_H2O <- remove_first_time(smeg3_df_bac_H2O)
smeg3_df_bac_H2O_spread <- spread(smeg3_df_bac_H2O, key = time, value = OD)
smeg3_df_bac_H2O <- filter(smeg3_df_bac_H2O_spread, smeg3_df_bac_H2O_spread$"8" < 1)

#make dataframe with only the antibiot ic wells
smeg3_df_pos_H2O <- filter(smeg3_df_H2O, well %in% H2O_positive_control_list)

#make inputs for mari3#####
mari3_df <- filter(df_metadata, bug == "marinum3")
mari3_df_DMSO <- filter(mari3_df , solvent == "DMSO")
mari3_df_H2O <- filter(mari3_df , solvent == "H2O")
mari3_df_MeOH <- filter(mari3_df , solvent == "MeOH")

mari3_df_library_DMSO <- filter(mari3_df_DMSO, well %notin% DMSO_filter_out)
mari3_df_library_DMSO <- remove_first_time(mari3_df_library_DMSO)
mari3_df_library_DMSO <- remove_first_time(mari3_df_library_DMSO)
mari3_df_library_DMSO_spread <- spread(mari3_df_library_DMSO, key = time, value = OD)
mari3_df_library_out <- filter(mari3_df_library_DMSO_spread, mari3_df_library_DMSO_spread$"8" > 0.35 | mari3_df_library_DMSO_spread$"1" > 0.15)
mari3_df_library_out <- mari3_df_library_out$well
mari3_df_library_DMSO <- filter(mari3_df_library_DMSO_spread, well %notin% mari3_df_library_out)
all_curves_solvent_new(mari3_df_library_DMSO)

#make dataframe with only the control wells
mari3_df_bac_DMSO <- filter(mari3_df_DMSO, well %in% DMSO_bac_control_list)
mari3_df_bac_DMSO <- remove_first_time(mari3_df_bac_DMSO)
mari3_df_bac_DMSO <- remove_first_time(mari3_df_bac_DMSO)
mari3_df_bac_DMSO_spread <- spread(mari3_df_bac_DMSO, key = time, value = OD)
mari3_df_bac_DMSO <- filter(mari3_df_bac_DMSO_spread, mari3_df_bac_DMSO_spread$"8" < 0.9 & mari3_df_bac_DMSO_spread$"8" > 0.15 & mari3_df_bac_DMSO_spread$"1" < 0.1)

#make dataframe with only the antibiot ic wells
mari3_df_pos_DMSO <- filter(mari3_df_DMSO, well %in% DMSO_positive_control_list)
mari_DMSO <- analysis_test_DMSO_2_inp(mari3_df_bac_DMSO, mari3_df_library_DMSO)
hist(mari_DMSO$cat$category)

#make similar dataframe but for MeOH
mari3_df_library_MeOH <- filter(mari3_df_MeOH, well %notin% MeOH_filter_out)
mari3_df_library_MeOH <- remove_first_time(mari3_df_library_MeOH)
mari3_df_library_MeOH_spread <- spread(mari3_df_library_MeOH, key = time, value = OD)
mari3_df_library_out <- filter(mari3_df_library_MeOH_spread, mari3_df_library_MeOH_spread$"8" > 1.25 |mari3_df_library_MeOH_spread$"1" > 0.25)
mari3_df_library_out <- mari3_df_library_out$well
mari3_df_library_MeOH <- filter(mari3_df_library_MeOH_spread, well %notin% mari3_df_library_out)
all_curves_solvent_new(mari3_df_library_MeOH)

#make dataframe with only the control wells
mari3_df_bac_MeOH <- filter(mari3_df_MeOH, well %in% MeOH_bac_control_list)
mari3_df_bac_MeOH <- remove_first_time(mari3_df_bac_MeOH)
mari3_df_bac_MeOH <- remove_first_time(mari3_df_bac_MeOH)
mari3_df_bac_MeOH_spread <- spread(mari3_df_bac_MeOH, key = time, value = OD)
mari3_df_bac_MeOH <- filter(mari3_df_bac_MeOH_spread, mari3_df_bac_MeOH_spread$"8" < 1 & mari3_df_bac_MeOH_spread$"1" < 0.25)
all_curves_solvent_new(mari3_df_bac_MeOH)
#make dataframe with only the antibiot ic wells
mari3_df_pos_MeOH <- filter(mari3_df_MeOH, well %in% MeOH_positive_control_list)

#make for H2O
mari3_df_library_H2O <- filter(mari3_df_H2O, well %notin% H2O_filter_out)
mari3_df_library_H2O <- remove_first_time(mari3_df_library_H2O)
mari3_df_library_H2O_spread <- spread(mari3_df_library_H2O, key = time, value = OD)
mari3_df_library_out <- filter(mari3_df_library_H2O_spread, mari3_df_library_H2O_spread$"8" > 1.25 | mari3_df_library_H2O_spread$"1" > 0.3)
mari3_df_library_out <- mari3_df_library_out$well
mari3_df_library_H2O <- filter(mari3_df_library_H2O_spread, well %notin% mari3_df_library_out)

#make dataframe with only the control wells
mari3_df_bac_H2O <- filter(mari3_df_H2O, well %in% H2O_bac_control_list)
mari3_df_bac_H2O <- remove_first_time(mari3_df_bac_H2O)
mari3_df_bac_H2O_spread <- spread(mari3_df_bac_H2O, key = time, value = OD)
mari3_df_bac_H2O <- filter(mari3_df_bac_H2O_spread, mari3_df_bac_H2O_spread$"8" < 0.35 & mari3_df_bac_H2O_spread$"1" < 0.1)
all_curves_solvent_new(mari3_df_bac_H2O)

#make dataframe with only the antibiot ic wells
mari3_df_pos_H2O <- filter(mari3_df_H2O, well %in% H2O_positive_control_list)

#make inputs for turi2#####
turi2_df <- filter(df_metadata, bug == "turicella2")
turi2_df_DMSO <- filter(turi2_df , solvent == "DMSO")
turi2_df_H2O <- filter(turi2_df , solvent == "H2O")
turi2_df_MeOH <- filter(turi2_df , solvent == "MeOH")

turi2_df_library_DMSO <- filter(turi2_df_DMSO, well %notin% DMSO_filter_out)
turi2_df_library_DMSO <- remove_first_time(turi2_df_library_DMSO)
turi2_df_library_DMSO <- remove_first_time(turi2_df_library_DMSO)
turi2_df_library_DMSO_spread <- spread(turi2_df_library_DMSO, key = time, value = OD)
turi2_df_library_out <- filter(turi2_df_library_DMSO_spread, turi2_df_library_DMSO_spread$"8" > 1.4 )
turi2_df_library_out <- turi2_df_library_out$well
turi2_df_library_DMSO <- filter(turi2_df_library_DMSO_spread, well %notin% turi2_df_library_out)
all_curves_solvent_new(turi2_df_library_DMSO)

#make dataframe with only the control wells
turi2_df_bac_DMSO <- filter(turi2_df_DMSO, well %in% DMSO_bac_control_list)
turi2_df_bac_DMSO <- remove_first_time(turi2_df_bac_DMSO)
turi2_df_bac_DMSO <- remove_first_time(turi2_df_bac_DMSO)
turi2_df_bac_DMSO_spread <- spread(turi2_df_bac_DMSO, key = time, value = OD)
turi2_df_bac_DMSO <- filter(turi2_df_bac_DMSO_spread, turi2_df_bac_DMSO_spread$"8" < 1.2)
all_curves_solvent_new(turi2_df_bac_DMSO)

#make dataframe with only the antibiot ic wells
turi2_df_pos_DMSO <- filter(turi2_df_DMSO, well %in% DMSO_positive_control_list)

#make similar dataframe but for MeOH
turi2_df_library_MeOH <- filter(turi2_df_MeOH, well %notin% MeOH_filter_out)
turi2_df_library_MeOH <- remove_first_time(turi2_df_library_MeOH)
turi2_df_library_MeOH <- remove_first_time(turi2_df_library_MeOH)
turi2_df_library_MeOH_spread <- spread(turi2_df_library_MeOH, key = time, value = OD)
turi2_df_library_out <- filter(turi2_df_library_MeOH_spread, turi2_df_library_MeOH_spread$"8" > 2 | turi2_df_library_MeOH_spread$"1" > 0.5)
turi2_df_library_out <- turi2_df_library_out$well
turi2_df_library_MeOH <- filter(turi2_df_library_MeOH_spread, well %notin% turi2_df_library_out)
all_curves_solvent_new(turi2_df_library_MeOH)

#make dataframe with only the control wells
turi2_df_bac_MeOH <- filter(turi2_df_MeOH, well %in% MeOH_bac_control_list)
turi2_df_bac_MeOH <- remove_first_time(turi2_df_bac_MeOH)
turi2_df_bac_MeOH <- remove_first_time(turi2_df_bac_MeOH)
turi2_df_bac_MeOH_spread <- spread(turi2_df_bac_MeOH, key = time, value = OD)
turi2_df_bac_MeOH <- filter(turi2_df_bac_MeOH_spread, turi2_df_bac_MeOH_spread$"8" > 0.5 & turi2_df_bac_MeOH_spread$"8" < 1.7)
all_curves_solvent_new(turi2_df_bac_MeOH_spread)
#make dataframe with only the antibiot ic wells
turi2_df_pos_MeOH <- filter(turi2_df_MeOH, well %in% MeOH_positive_control_list)

#make for H2O
turi2_df_library_H2O <- filter(turi2_df_H2O, well %notin% H2O_filter_out)
turi2_df_library_H2O <- remove_first_time(turi2_df_library_H2O)
turi2_df_library_H2O <- remove_first_time(turi2_df_library_H2O)
turi2_df_library_H2O_spread <- spread(turi2_df_library_H2O, key = time, value = OD)
turi2_df_library_out <- filter(turi2_df_library_H2O_spread, turi2_df_library_H2O_spread$"8" > 1.5 | turi2_df_library_H2O_spread$"1" > 0.5)
turi2_df_library_out <- turi2_df_library_out$well
turi2_df_library_H2O <- filter(turi2_df_library_H2O_spread, well %notin% turi2_df_library_out)
all_curves_solvent_new(turi2_df_library_H2O)

#make dataframe with only the control wells
turi2_df_bac_H2O <- filter(turi2_df_H2O, well %in% H2O_bac_control_list)
turi2_df_bac_H2O <- remove_first_time(turi2_df_bac_H2O)
turi2_df_bac_H2O <- remove_first_time(turi2_df_bac_H2O)
turi2_df_bac_H2O_spread <- spread(turi2_df_bac_H2O, key = time, value = OD)
turi2_df_bac_H2O <- filter(turi2_df_bac_H2O_spread, turi2_df_bac_H2O_spread$"8" < 1.25 & turi2_df_bac_H2O_spread$"8" > 0.5)
all_curves_solvent_new(turi2_df_bac_H2O)

#make dataframe with only the antibiot ic wells
turi2_df_pos_H2O <- filter(turi2_df_H2O, well %in% H2O_positive_control_list)

turi2_DMSO <- analysis_test_DMSO_2_inp(turi2_df_bac_DMSO, turi2_df_library_DMSO)
turi2_MeOH <- analysis_test_DMSO_2_inp(turi2_df_bac_MeOH, turi2_df_library_MeOH)
turi2_H2O <- analysis_test_DMSO_2_inp(turi2_df_bac_H2O, turi2_df_library_H2O)


###FILTER OUT TIME POINT 1####
#input = whatever dataframe you want to remove the first time point from
remove_first_time <- function(input){
  input_1 <- filter(input, time != "1")
  input_1$time <- as.numeric(input_1$time)
  input_1$time <- input_1$time-1 
  return(input_1)
}

remove_9_time <- function(input){
  input <- filter(input, time != "9")
  return(input)
}



analysis_test_DMSO_2_inp <- function(bac_df, lib_df){
  #format in usable shape for library
  arth_df_library_DMSO <- lib_df
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
  arth_df_bac_DMSO <- bac_df
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
  
  
  #same but for pos controls
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
  results <- list(bac_od = arth_df_bac_DMSO_od_vector, bac_spline = arth_df_bac_DMSO_spline_vector, cat = ttest_output_arth_DMSO) 
  return(results)
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


#same but for pos controls
arth_df_pos_DMSO <- spread(arth_df_pos_DMSO, key = time, value = OD)
#filter values above OD 4.5 - example that could work
#arth_df_pos_DMSO <- arth_df_pos_DMSO[arth_df_pos_DMSO$"8" < 4.5]
###################TEST FILTER###########################
arth_df_pos_DMSO$diff <- arth_df_pos_DMSO$"8" - arth_df_pos_DMSO$"1"
arth_df_pos_DMSO_od_vector <- arth_df_pos_DMSO$diff

#vector for spline
arth_df_pos_DMSO_splines <- select(arth_df_pos_DMSO, c(6:13))
arth_pos_time_DMSO <- c(0:7)
arth_pos_lm_DMSO <- apply(arth_df_pos_DMSO_splines, 1, function(x) lm(x ~ lspline(arth_pos_time_DMSO, knots = knots)))
arth_df_pos_DMSO_splines_coeff <- as.data.frame(unlist(lapply(arth_pos_lm_DMSO, function(x) coef(x)[2])))
#combine
arth_df_pos_DMSO$spline <- arth_df_pos_DMSO_splines_coeff
arth_df_pos_DMSO_spline_vector <- arth_df_pos_DMSO$spline

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
results <- list(bac_od = arth_df_bac_DMSO_od_vector, bac_spline = arth_df_bac_DMSO_spline_vector, pos_spline = arth_df_pos_DMSO_spline_vector, pos_od = arth_df_pos_DMSO_od_vector, cat = ttest_output_arth_DMSO) 
return(results)
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
  
  #same but for pos controls
  arth_df_pos_H2O <- spread(arth_df_pos_H2O, key = time, value = OD)
  #filter values above OD 4.5 - example that could work
  #arth_df_pos_H2O <- arth_df_pos_H2O[arth_df_pos_H2O$"8" < 4.5]
  ###################TEST FILTER###########################
  arth_df_pos_H2O$diff <- arth_df_pos_H2O$"8" - arth_df_pos_H2O$"1"
  arth_df_pos_H2O_od_vector <- arth_df_pos_H2O$diff
  
  #vector for spline
  arth_df_pos_H2O_splines <- select(arth_df_pos_H2O, c(6:13))
  arth_pos_time_H2O <- c(0:7)
  arth_pos_lm_H2O <- apply(arth_df_pos_H2O_splines, 1, function(x) lm(x ~ lspline(arth_pos_time_H2O, knots = knots)))
  arth_df_pos_H2O_splines_coeff <- as.data.frame(unlist(lapply(arth_pos_lm_H2O, function(x) coef(x)[2])))
  #combine
  arth_df_pos_H2O$spline <- arth_df_pos_H2O_splines_coeff
  arth_df_pos_H2O_spline_vector <- arth_df_pos_H2O$spline
  
  
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
  results <- list(bac_od = arth_df_bac_H2O_od_vector, bac_spline = arth_df_bac_H2O_spline_vector, pos_spline = arth_df_pos_H2O_spline_vector, pos_od = arth_df_pos_H2O_od_vector, cat = ttest_output_arth_H2O) 
  return(results)
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
  
  #same but for positive controls
  arth_df_pos_MeOH <- spread(arth_df_pos_MeOH, key = time, value = OD)
  #filter values above OD 4.5 - example that could work
  #arth_df_pos_MeOH <- arth_df_pos_MeOH[arth_df_pos_MeOH$"8" < 4.5]
  ###################TEST FILTER###########################
  arth_df_pos_MeOH$diff <- arth_df_pos_MeOH$"8" - arth_df_pos_MeOH$"1"
  arth_df_pos_MeOH_od_vector <- arth_df_pos_MeOH$diff
  
  #vector for spline
  arth_df_pos_MeOH_splines <- select(arth_df_pos_MeOH, c(6:13))
  arth_pos_time_MeOH <- c(0:7)
  arth_pos_lm_MeOH <- apply(arth_df_pos_MeOH_splines, 1, function(x) lm(x ~ lspline(arth_pos_time_MeOH, knots = knots)))
  arth_df_pos_MeOH_splines_coeff <- as.data.frame(unlist(lapply(arth_pos_lm_MeOH, function(x) coef(x)[2])))
  #combine
  arth_df_pos_MeOH$spline <- arth_df_pos_MeOH_splines_coeff
  arth_df_pos_MeOH_spline_vector <- arth_df_pos_MeOH$spline
  
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
  ttest_output_arth_MeOH <- merge(arth_MeOH_library, ttest_output_arth_MeOH, by = "well" )
  results <- list(bac_od = arth_df_bac_MeOH_od_vector, bac_spline = arth_df_bac_MeOH_spline_vector, pos_spline = arth_df_pos_MeOH_spline_vector, pos_od = arth_df_pos_MeOH_od_vector, cat = ttest_output_arth_MeOH) 
  return(results)
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




#Z' calculation
z_test <- function(x){
  df_bac_DMSO <- filter(x, well %in% DMSO_bac_control_list)
  st_d_bac <- filter(df_bac_DMSO, time == "8")
  st_d_bac_3 <- 3 * sd(st_d_bac$OD)
  mean_bac <- mean(st_d_bac$OD)
  
  df_pos_DMSO <- filter(x, well %in% DMSO_positive_control_list)
  st_d_pos <- filter(df_pos_DMSO, time == "8")
  st_d_pos_3 <- 3 * sd(st_d_pos$OD)
  mean_pos <- mean(st_d_pos$OD)
  
  
  z <- 1 - ((st_d_pos_3 + st_d_bac_3)/(abs(mean_pos - mean_bac)))
  return(z)
}


s <- (mean_pos - st_d_pos_3) - (mean_bac + st_d_bac_3) 
r <- mean_pos - mean_bac
s/r

ggplot(st_d_bac, aes(x=OD)) + geom_density()
ggplot(st_d_pos, aes(x=OD)) + geom_density()



df_bac_DMSO <- filter(arth_df_DMSO, well %in% DMSO_bac_control_list)
df_bac_DMSO <- spread(df_bac_DMSO, key = time, value = OD)
  df_bac_DMSO <- filter(df_bac_DMSO, df_bac_DMSO$"8" < 3)
st_d_bac <- df_bac_DMSO$"8"
st_d_bac_3 <- 3 * sd(st_d_bac)
mean_bac <- mean(st_d_bac)

df_pos_DMSO <- filter(arth_df_DMSO, well %in% DMSO_positive_control_list)
df_pos_DMSO <- spread(df_pos_DMSO, key = time, value = OD)
df_pos_DMSO <- filter(df_pos_DMSO, df_pos_DMSO$"8" < 0.5)
st_d_pos <- df_pos_DMSO$"8"
st_d_pos_3 <- 3 * sd(st_d_pos)
mean_pos <- mean(st_d_pos)


z <- 1 - ((st_d_pos_3 + st_d_bac_3)/(abs(mean_pos - mean_bac)))

#dOD600 and spline
z_d600_spline <- function(x){

sd_bac_d600 <- sd(x$bac_od)
mean_bac_d600 <- mean(x$bac_od)
sd_pos_d600 <- sd(x$pos_od)
mean_pos_d600 <- mean(x$pos_od)

sd_bac_spline <- sd(x$bac_spline)
mean_bac_spline <- mean(x$bac_spline)
sd_pos_spline <- sd(x$pos_spline)
mean_pos_spline <- mean(x$pos_spline)

z_d600 <- 1 - ((3*sd_pos_d600 + 3*sd_bac_d600)/abs(mean_pos_d600 - mean_bac_d600))
z_spline <- 1 - ((3*sd_pos_spline + 3*sd_bac_spline)/abs(mean_pos_spline - mean_bac_spline))
results <- list(z_d600 = z_d600, z_spline = z_spline)
return(results)
}

#functions to run for z-prime (unfiltered)
rhodo_DMSO <-analysis_test_DMSO(rhodo_df_DMSO)
rhodo_DMSO_z <- z_d600_spline(rhodo_DMSO)

rhodo_MeOH <-analysis_test_MeOH(rhodo_df_MeOH)
rhodo_MeOH_z <- z_d600_spline(rhodo_MeOH)

rhodo_H2O <- analysis_test_H2O(rhodo_df_H2O)
rhodo_H2O_z <- z_d600_spline(rhodo_H2O)

mari3_DMSO <-analysis_test_DMSO(mari3_df_DMSO)
mari3_DMSO_z <- z_d600_spline(mari3_DMSO)

mari3_MeOH <-analysis_test_MeOH(mari3_df_MeOH)
mari3_MeOH_z <- z_d600_spline(mari3_MeOH)

mari3_H2O <- analysis_test_H2O(mari3_df_H2O)
mari3_H2O_z <- z_d600_spline(mari3_H2O)

turi2_MeOH_out <- turi2_MeOH$cat
colnames(turi2_MeOH_out) <- c("Well", "cat")
turi2_MeOH_out_test <- merge(platemap_MeOH, turi2_MeOH_out, by = "Well" )
write.csv(turi2_MeOH_out_test, "turi2_MeOH_categories.csv")

turi2_DMSO_out <- turi2_DMSO$cat
colnames(turi2_DMSO_out) <- c("Well", "cat")
turi2_DMSO_out_test <- merge(platemap_DMSO, turi2_DMSO_out, by = "Well" )
write.csv(turi2_DMSO_out_test, "turi2_DMSO_categories.csv")

turi2_H2O_out <- turi2_H2O$cat
colnames(turi2_H2O_out) <- c("Well", "cat")
turi2_H2O_out_test <- merge(platemap_H2O, turi2_H2O_out, by = "Well" )
write.csv(turi2_H2O_out_test, "turi2_H2O_categories.csv")
