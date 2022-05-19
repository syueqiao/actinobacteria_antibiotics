#As of May 2022 this is the working version

#load packages that are needed
library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(lspline)

#set wd
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/strains_data")
#set this useful function
`%notin%` <- Negate(`%in%`)

#import platemap
#see supporting files for this document
platemap <- read.csv("../plate_maps.csv", header = TRUE)

platemap_DMSO <- select(platemap, Well, DMSO)
DMSO_media_list <- platemap_DMSO[platemap_DMSO$DMSO %like% "Empty",]  
DMSO_media_list <- as.vector(DMSO_media_list$Well)

DMSO_bac_control_list <- platemap_DMSO[platemap_DMSO$DMSO %like% "Bacteria only", ]  
DMSO_bac_control_list <- as.vector(DMSO_bac_control_list$Well)

DMSO_positive_control_list <- platemap_DMSO[platemap_DMSO$DMSO %like% "Positive", ]
DMSO_positive_control_list <- as.vector(DMSO_positive_control_list$Well)


#set up loop for importing the metadata based ON THE DIRECTORY NAMES!!!
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

#create list of controls and stuff to filter out
DMSO_filter_out <- c(DMSO_media_list, DMSO_positive_control_list, DMSO_bac_control_list)

arth_df <- filter(df_metadata, bug == "ArthBac")
arth_df_DMSO <- filter(arth_df, solvent == "DMSO")

#make dataframe with only the library wells
arth_df_library_DMSO <- filter(arth_df_DMSO, well %notin% DMSO_filter_out)
#make dataframe with only the control wells
arth_df_bac_DMSO <- filter(arth_df_DMSO, well %in% DMSO_bac_control_list)
#make dataframe with only the antibiotic wells
arth_df_pos_DMSO <- filter(arth_df_DMSO, well %in% DMSO_positive_control_list)

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
ttest_output_arth_DMSO <- merge(arth_DMSO_library, ttest_output_arth_DMSO, by = "well" )

#export data
path <- 'C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics/category_outputs'
write.csv(ttest_output_arth_DMSO, file.path(path, "arth_DMSO_categories.csv"), row.names = F)

#plot prelim figure
ggplot(ttest_output_arth_DMSO, aes(x=category, fill=..x..)) +
  geom_bar() +
  xlab("Category") +
  ylab("Count") +
  ggtitle("Effect on Arthrobacter growth of DMSO-solvent compounds") +
  scale_x_continuous(breaks = 1:9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = '#5AAA46', high = '#825CA6', breaks = c(1:9), labels = c(1:9), name = "Categories")



##testing the antibiotic only control stuff
###################################################################################
arth_df_pos_DMSO <- spread(arth_df_pos_DMSO, key = time, value = OD)
#find diff from t8-t1
arth_df_pos_DMSO$diff <- arth_df_pos_DMSO$"8" - arth_df_pos_DMSO$"1"


#find spline from first two time points
#select the values only
arth_df_pos_DMSO_splines <- select(arth_df_pos_DMSO, c(6:13))

#make time vector of right length
arth_pos_time_DMSO <- c(0:7)

#find splines and extract coefficient
knots <- c(3)
arth_lm_DMSO <- apply(arth_df_pos_DMSO_splines, 1, function(x) lm(x ~ lspline(arth_pos_time_DMSO, knots = knots)))
arth_df_pos_DMSO_splines_coeff <- as.data.frame(unlist(lapply(arth_lm_DMSO, function(x) coef(x)[2])))
colnames(arth_df_pos_DMSO_splines_coeff) <- c("spline")

#combine with dig dataframe
arth_df_pos_DMSO$spline <- arth_df_pos_DMSO_splines_coeff

#spread the data for t tests
arth_df_pos_DMSO_spread_spline <- spread(arth_df_pos_DMSO, key = well, value = spline)
pos_cols1 <- ncol(arth_df_pos_DMSO)-1
pos_cols2 <- ncol(arth_df_pos_DMSO_spread_spline)
arth_df_pos_DMSO_spread_spline <- select(arth_df_pos_DMSO_spread_spline, c(pos_cols1:pos_cols2))
arth_df_pos_DMSO_spread_spline <- data.frame(lapply(arth_df_pos_DMSO_spread_spline, na.omit))

arth_df_pos_DMSO_spread_diff <- spread(arth_df_pos_DMSO, key = well, value = diff)
arth_df_pos_DMSO_spread_diff <- select(arth_df_pos_DMSO_spread_diff, c(pos_cols1:pos_cols2))
arth_df_pos_DMSO_spread_diff <- data.frame(lapply(arth_df_pos_DMSO_spread_diff, na.omit))

#get vectors for the bac only 



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
arth_ttests_od_DMSO <- lapply(arth_df_pos_DMSO_spread_diff, function(x) wilcox.test(x, arth_df_bac_DMSO_od_vector, exact = T)$p.value)
arth_ttests_od_DMSO <- as.data.frame(unlist(arth_ttests_od_DMSO))
colnames(arth_ttests_od_DMSO) <- c("OD diff significance")
#wilcox.test(arth_df_bac_od_vector, arth_df_time_spread_diff$C15, exact = T)

#this is for the splines
#wilcox.test(arth_df_bac_splines_vector, arth_df_time_spread_spline$C15, exact = T)
arth_ttests_spline_DMSO <- lapply(arth_df_pos_DMSO_spread_spline, function(x) wilcox.test(x, arth_df_bac_DMSO_spline_vector, exact = T)$p.value)
arth_ttests_spline_DMSO <- as.data.frame(unlist(arth_ttests_spline_DMSO))
colnames(arth_ttests_spline_DMSO) <- c("spline diff significance")

#combine the two together
arth_ttests_combined_DMSO <- cbind(arth_ttests_od_DMSO, arth_ttests_spline_DMSO)

#rownames to column
arth_ttests_combined_DMSO <- rownames_to_column(arth_ttests_combined_DMSO, "well")

#compare?
#make it so that you can compare values
arth_ttests_combined_plus_data_DMSO <- merge(arth_ttests_combined_DMSO, arth_df_pos_DMSO, by = "well")
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

#QC and visualization of the outputs
hist(arth_ttests_combined_plus_data_DMSO$category, breaks = 12)
table(arth_ttests_combined_plus_data_DMSO$category)
arth_ttests_histo <- c(arth_ttests_combined_plus_data_DMSO$category)
ggplot(ttest_output_arth_DMSO, aes(x=category, fill=..x..)) +
  geom_bar() +
  xlab("Category") +
  ylab("Count") +
  ggtitle("Effect on Arthrobacter growth of DMSO-solvent compounds") +
  scale_x_continuous(breaks = 1:9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = '#5AAA46', high = '#825CA6', breaks = c(1:9), labels = c(1:9), name = "Categories")

#make nicer on the eyes
ttest_output_arth_DMSO <- select(arth_ttests_combined_plus_data_DMSO, well, category)
ttest_output_arth_DMSO <- ttest_output_arth_DMSO[!duplicated(ttest_output_arth_DMSO), ]

arth_DMSO_pos <- platemap_DMSO
names(arth_DMSO_pos)[names(arth_DMSO_pos) == 'Well'] <- 'well'

ttest_output_arth_DMSO <- merge(arth_DMSO_pos, ttest_output_arth_DMSO, by = "well" )
#write.csv(ttest_output_arth_DMSO, "arth_DMSO_categories.csv")


