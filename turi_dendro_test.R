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

#read the file test1.smiles <- you'll need to set up anything input into the right format
turi_all_smiles <- read_tsv("turi_all_smiles.txt")
turi_all_smiles <- turi_all_smiles[!duplicated(turi_all_smiles$Chemical),]

turi_smiles_NA <- turi_all_smiles_test[is.na(turi_all_smiles_test$SMILE),]
turi_smiles_NA <- left_join(turi_smiles_NA, turi_all_NA, by = "Chemical")
turi_smiles_NA <- subset(turi_smiles_NA, select = -c(SMILE.x))
write_tsv(turi_smiles_NA, "turi_smiles_NA_new.tsv")


turi_all_NA <- read_tsv("turi_NAs.txt")
turi_all_NA <- turi_all_NA[!duplicated(turi_all_NA$Chemical),]
turi_all_smiles <- na.omit(turi_all_smiles)

turi_smiles <- rbind(turi_all_smiles, turi_all_NA)
turi_smiles <- rev(turi_smiles)
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

hc <- hclust(as.dist(distmat), method="ward") 
hc[["labels"]] <- cid(apset) # Assign correct item labels 
plot(as.dendrogram(hc), edgePar=list(col=50, lwd=2), horiz=T)

png(filename = "testy5.png", res = 300,
    width = 10000, height = 10000)
t <- plot(as.phylo(hc), type = "unrooted", cex = 0.2, lab4ut = "axial", 
     no.margin = TRUE)
plot(hc, hang = -1, cex = 0.3)
dev.off()

t <- as.phylo(hc)
write.tree(t, "test_tree2.txt")

p <- plot(as.phylo(hc), type = "unrooted", cex = 0.2, lab4ut = "axial", label.offset = 1,
     no.margin = TRUE)

png(filename = "testy2.png", res = 300,
    width = 2000, height = 2000)

plot(as.phylo(hc), type = "unrooted", cex = 0.4,
     no.margin = TRUE)
dev.off()

#use the hclust function in R to cluster them based on distance matrix
sweet_dendro <- as.dendrogram(hclust(d = dist(x = sweet_matrix)))
#plot w/ ggdendrogram
dendro_plot <- ggdendrogram(data = sweet_dendro, rotate = TRUE)
dendro_plot
#import the categories file. I made a few adjusments to this and you might have to as well
turi_Me_cat <- read.csv("turi_Me_categories.csv")
turi_Me_cat <- select(turi_Me_cat, c("index", "well", "category", "bac"))


#order the y-axis in the order of the dendrogram tree   
sweet_order <- order.dendrogram(sweet_dendro)
turi_Me_cat$index <- factor(x = turi_Me_cat$index,
                            levels = turi_Me_cat$index[sweet_order], 
                            ordered = TRUE)



#repeat lines 56-58 for each bac and solvent

#plot heatmap of categories
heatmap_plot <- ggplot(data = turi_Me_cat, aes(x = bac, y = index)) +
  geom_tile(aes(fill = category)) +
  scale_fill_distiller(palette = "PRGn") +
  theme(axis.text.y = element_text(size = 6), legend.position = "top")
#show
heatmap_plot    

#use grid function to slap the two figures together. need to manually adjust to make things line up
grid.newpage()
print(heatmap_plot, 
      vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro_plot, 
      vp = viewport(x = 0.90, y = 0.47, width = 0.2, height = 1.03))

#repeat lines 62-75 for each bac and solvent

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

colnames(labels_df) <- c("solv")
turi2_all_out_fuzzy_jw <- stringdist_join(labels_df, turi2_all_out_merge, method = "jw", by = "solv", max_dist = 99, distance_col = 'dist') %>% group_by(solv.x) %>% slice_min(order_by = dist, n=1)

write.csv(turi2_all_out_fuzzy_jw, "turi_out_fuzzy_jw.csv")

#DO NOT EDIT THE OUTPUT FILE HERE IT WILL BE SCREWED UP AND HARD TO UNDO
turi_all_out_fuzzy_jw_edited <- read.csv("turi_out_fuzzy_jw_edited.csv")

#tree with only inhibitors
turi2_inhib_out_merge <- filter(turi_all_out_fuzzy_jw_edited, category == "9" | category == "8")
inhibs <- data.frame(turi2_inhib_out_merge$solv.x)
colnames(inhibs) <- c("Chemical")
inhibs_smile <- stringdist_join(inhibs, turi_smiles, method = "jw", by = "Chemical", max_dist = 99, distance_col = "dist") %>% group_by(Chemical.x) %>% slice_min(order_by = dist, n=1)
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
