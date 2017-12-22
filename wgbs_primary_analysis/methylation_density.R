# remaking plot from Jensen et al 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0645-x#MOESM1
# using my implementation of ENCODE WGBS pipeline 
# 15 Nov 17 
# author <christa.caggiano@ucsf.edu> 

# the environment better be empty because this is a hella lot of data 
rm(list = ls())
setwd("~/Documents/UCSF_year1/Zaitlen-rotation1/trash/")

# data.table library and fread to read in large text files easier 
library(data.table)

# load placenta samples 
plac = fread("Placenta_1_merged.txt")
plac2 = fread("placenta2_merged.txt")
plac3 = fread("Placenta_4_merged.txt")

# load non-pregnant cfDNA samples 
nonpreg1 = fread("NonPregnant.ccfDNA_0_merged.txt")
nonpreg2 = fread("NonPregnant.ccfDNA_4_merged.txt")
nonpreg3 = fread("NonPregnant.ccfDNA_3_merged.txt")

n = subset(nonpreg3, nonpreg3$V1=="chr2" & nonpreg3$V2 > 79347400 & nonpreg3$V2 < 79347595)

# load buffy coat samples 
buffy1 = fread("Buffy.Coat_2_merged.txt")
buffy2 = fread("Buffy.Coat_3_merged.txt")
buffy3 = fread("Buffy.Coat_6_merged.txt")

# filter by coverage - 100 as a cut-off relatively arbitrary, but replicates 
# jensen paper well 
buffy2_filt = subset(buffy2, buffy2$V5 > 100 )
buffy1_filt = subset(buffy1, buffy1$V5 > 100)
buffy3_filt = subset(buffy3, buffy3$V5 > 100)

nonpreg1_filt = subset(nonpreg1, nonpreg1$V5 > 100)
nonpreg2_filt = subset(nonpreg2, nonpreg2$V5 >100)
nonpreg3_filt = subset(nonpreg3, nonpreg3$V5 > 100)

plac1_filt = subset(plac, plac$V5 > 100)
plac2_filt = subset(plac2, plac2$V5 > 100)
plac3_filt = subset(plac3, plac3$V5 > 100)

# merge all replicates of the samples together 
# just straight pooled together because we are looking @ averages 
buffy = rbind(buffy2_filt, buffy1_filt, buffy3_filt)
plac = rbind(plac1_filt, plac2_filt, plac3_filt)
nonpreg = rbind(nonpreg1_filt, nonpreg2_filt, nonpreg3_filt)

# remove edges - just testing this, trying to replicate plot 
buffy = subset(buffy, buffy$V6>0 & buffy$V6 <100)
plac = subset(plac, plac$V6>0 & plac$V6 <100)
nonpreg = subset(nonpreg, nonpreg$V6>0 & nonpreg$V6 <100)

# sets the number of bins for the histogram 
bins = seq(0, 100)

# produce histogram using R base graphics 
png(filename="jensenhist.png")
hist(plac$V6, col = "Red", freq=FALSE, breaks=bins, ylim=c(0, 0.08), xlab = "% methylation", ylab="proportion", main="CpG methlyation")
hist(buffy$V6, col = "Blue", breaks=bins, freq=FALSE, add=T)
hist(nonpreg$V6, col = "Green", breaks=bins, freq=FALSE, add=T)
dev.off()

# use ggplot to do a density curve plot 
library(ggplot2)

# stupid data manipulation 
vec1 = data.frame(x=plac$V6)
vec2 = data.frame(x=buffy$V6)
vec3 = data.frame(x=nonpreg$V6)

# preliminary density curve 
ggplot() + geom_density(aes(x=x), colour="red", data=vec1) + 
  geom_density(aes(x=x), colour="blue", data=vec2) +  geom_density(aes(x=x), colour="green", data=vec3)  

# clear all the junk out of the memory just in case 
rm(list = ls())

