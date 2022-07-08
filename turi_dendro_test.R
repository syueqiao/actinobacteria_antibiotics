#script to make dendrograms for the turi data
setwd("C:/Users/Jessica Shen/Desktop/actinobacteria_antibiotics")
#might need to install a few of the packages. every package should have google-able installation instructions

library("ChemmineR")
library("ChemmineOB")
library("tidyverse")
library("ggplot2")
library("ggdendro")
library("grid")
library("RColorBrewer")
library("dendextend")
library("ape")
library("fuzzyjoin")

#read the file test1.smiles <- you'll need to set up anything input into the right format
turi_all_smiles <- read_tsv("turi_all_smiles.txt")

#remove duplicates because pubchem search returns duplicates
turi_all_smiles <- turi_all_smiles[!duplicated(turi_all_smiles$Chemical),]

#subset out the chemicals pubchem was not able to find smiles for
turi_smiles_NA <- turi_all_smiles_test[is.na(turi_all_smiles_test$SMILE),]
write_tsv(turi_smiles_NA, "turi_smiles_NA_new.tsv")

#manually edit in the smiles that pubchem was not able to get
turi_all_NA <- read_tsv("turi_NAs.txt")
turi_all_NA <- turi_all_NA[!duplicated(turi_all_NA$Chemical),]
turi_all_smiles <- na.omit(turi_all_smiles)

#combine two dataframes together 
turi_smiles <- rbind(turi_all_smiles, turi_all_NA)
#reverse dataframes since chemmineR assumes the first column is smiles, and the second column is the identifier
turi_smiles <- rev(turi_smiles)
#write out to import
write_tsv(turi_smiles, "turi_smiles_import.txt")
#need to manually add some files that were not able to be automatically retrieved######

#assume Ergotheoneine is Ergothioneine
#assume lanthium is Lanthanum (III) carbonate hydrate
#assume Norfloxicin is Norfloxacin
#DL-Fluorocitrate barium salt is DL-Fluorocitric barium salt
#some other typos, weird

smiset <- read.SMIset("turi_smiles_import.txt")
#convert to chemminer working format
sdfset <- smiles2sdf(smiset)
#convert to apset format for clustering
apset <- sdf2ap(sdfset)
#important!!! run this to make sure that the number of compounds is what you expect, eg, dmso has 142 compounds
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
write.tree(t, "test_tree2.txt")

#do stuff that allows a category to be assignd to each label, basically, 

platemap_DMSO_merge <- platemap_DMSO
colnames(platemap_DMSO_merge) <- c("well", "solv")
turi2_DMSO_out_merge <- merge(turi2_DMSO$cat, platemap_DMSO_merge, by = "well")

platemap_H2O_merge <- platemap_H2O
colnames(platemap_H2O_merge) <- c("well", "solv")
turi2_H2O_out_merge <- merge(turi2_H2O$cat, platemap_H2O_merge, by = "well")


platemap_MeOH_merge <- platemap_MeOH
colnames(platemap_MeOH_merge) <- c("well", "solv")
turi2_MeOH_out_merge <- merge(turi2_MeOH$cat, platemap_MeOH_merge, by = "well")

turi2_all_out_merge <- rbind(turi2_MeOH_out_merge, turi2_H2O_out_merge, turi2_DMSO_out_merge)

#perform "fuzzy merge"
#make dataframe with the labels
labels_df <- data.frame(cid(apset))
colnames(labels_df) <- c("solv")
turi2_all_out_fuzzy_jw <- stringdist_join(labels_df, turi2_all_out_merge, method = "jw", by = "solv", max_dist = 99, distance_col = 'dist') %>% group_by(solv.x) %>% slice_min(order_by = dist, n=1)

write.csv(turi2_all_out_fuzzy_jw, "turi_out_fuzzy_jw.csv")

#DO NOT EDIT THE OUTPUT FILE HERE IT WILL BE SCREWED UP AND HARD TO UNDO
#basically, go through manually and make sure that the fuzzy merge was done correctly
#sort by distance, the higher the "distance" the more likely that the matching was done incorrectly

#import edited version
turi_all_out_fuzzy_jw_edited <- read.csv("turi_out_fuzzy_jw_edited.csv")


#from this point on you can filter and make trees for each category assigned
#tree with only inhibitors
turi2_inhib_out_merge <- filter(turi_all_out_fuzzy_jw_edited, category == "9" | category == "8")

#grab chemical labels that are the inhibitors
inhibs <- data.frame(turi2_inhib_out_merge$solv.x)
colnames(inhibs) <- c("Chemical")

#once again, fuzzy merge with the list of smiles, so that you can get the smiles for each chemical
#i fuzzy merge here because the cid(apset) function seems to remove the spaces
#this is generally very consistent and won't need to manually change much
inhibs_smile <- stringdist_join(inhibs, turi_smiles, method = "jw", by = "Chemical", max_dist = 99, distance_col = "dist") %>% group_by(Chemical.x) %>% slice_min(order_by = dist, n=1)
#clean up format for putting into clustering program
inhibs_smile <- select(inhibs_smile, Chemical.x, SMILE)
inhibs_smile <- rev(inhibs_smile)
write_tsv(inhibs_smile, "turi_inhibs_smiles_import.txt")

#need to delete header manually
smiset_inhib <- read.SMIset("turi_inhibs_smiles_import.txt")
#convert to chemminer working format
sdfset_inhib <- smiles2sdf(smiset_inhib)
#convert to apset format for clustering
apset_inhib <- sdf2ap(sdfset_inhib)
#important!!! run this to make sure that the number of compounds is what you expect, eg, dmso has 142 compounds
apset_inhib

#this is the clustering command from ChemmineR. I'm not exactly sure what the cutoffs do
#maybe take a look?
#this will save the file "distmat.rda" into your working directory
clusters_inhib <- cmp.cluster(db=apset_inhib, cutoff = c(0.7, 0.8, 0.9),
                        save.distances="distmat_inhib.rda")

#put distamat.rda into R environment
load("distmat_inhib.rda") 

#convert to matrix format
sweet_matrix_inhib <- as.matrix(distmat)

hc_inhib <- hclust(as.dist(distmat), method="ward") 
hc_inhib[["labels"]] <- cid(apset_inhib) # Assign correct item labels 

t_inhib <- as.phylo(hc_inhib)
write.tree(t_inhib, "test_tree_inhib.txt")
