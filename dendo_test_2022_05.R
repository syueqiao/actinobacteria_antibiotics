setwd("C:/Users/Lauren/Desktop/R")
#might need to install a few of the packages. every package should have google-able installation instructions
library("ChemmineR")
library("ChemmineOB")
library("tidyverse")
library("ggplot2")
library("ggdendro")
library("grid")
library("RColorBrewer")

#read the file test1.smiles <- you'll need to set up anything input into the right format
smiset <- read.SMIset("Me_smiles.txt")
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