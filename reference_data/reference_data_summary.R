# script to summarize results of reference methylation datasets
# dec 2017 
# author <christa.caggiano@ucsf.edu> 

############################# environment #############################################

rm(list=ls()) # removes anything in memory 

# requires "wgbs_secondary_data_key.txt" and "pnas_sites.txt" to be in wd 
# sets working directory 
# setwd("~/Desktop/zaitlen_lab_desktop/")

library(wesanderson) # library of wes anderson themed color palettes, for funsies
require(gplots) # updates to R base graphics 

############# analysis of metadata associated with methylation datasets ################

wgbs = read.csv("wgbs_secondary_data_key.txt", sep="\t")

wgbs$Tissue=as.factor(wgbs$Tissue) # makes sure data loaded properly 
wgbs_no_blood = subset(wgbs, wgbs$Tissue!=" whole blood") # removes whole blood, too many are included to be balanced

# plots tissue types 
pal = wes_palette("Cavalcanti", nrow(table(s_no_blood$Tissue)), type = "continuous") # makes a wes anderson color palette for barplot
op <- par(mar = c(10,4,4,2) + 0.1) # sets the graphing parameters so that there is extra whitespace around the graph, good for long x-labels on the barplot
barplot(table(s_no_blood$Tissue), col=pal, las=2, cex.names=0.7, ylim=c(0, 30)) # plots (las rotates text)
title(main="compiled tissue types", ylab="counts", las=2) # adds title 

# creates a bar plot of gender breakdown in methylation datasets 
table(s_no_blood$Gender)
pie(as.factor(s_no_blood$Gender), labels=c("Male", "Female", "unknown"))

############# analysis reference data percent methylated at select pnas sites ################
pnas = read.table("pnas_sites.txt", header=T) # read data 

pnas = pnas[4:40] # subset to just % methylated columns 
pnas = na.omit(pnas) # omit rows with an NAs 

heatmap.2(cov(pnas), Rowv=NA, Colv=NA, col = greenred(256), trace="none", scale="column", dendrogram="none") # make a heatmap of covariance between sites




