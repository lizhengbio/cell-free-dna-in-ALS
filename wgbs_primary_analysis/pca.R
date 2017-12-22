# merging percent methylation files on CpG files
# pca analysis on als meth data 
# dec 4 2017
# author <christa.caggiano@ucsf.edu> 

############################# environment ##############################

rm(list = ls()) # clear memory

library(grDevices) # library for plotting 
library(ggplot2) # also for plotting 

# source("https://bioconductor.org/biocLite.R")
# biocLite("pcaMethods")
library(pcaMethods) # library for performing pca with missing data 


# set working directory
# setwd("~/Desktop/")

# merges chr and start text 
make.key <- function(x) {
  x$key <-paste(x$V1, x$V2)
  return(x)
}

# calculate percent methylation 
percent <- function(x) { 
  x$percent <- x$V3 / x$V5 * 100
  return(x)
}

# keep specific columns to merge with 
keep <- function(x) {
  x <- x[, c("percent", "key")]
  return(x)
}

# subset to methylation sites greater than n
filter <- function(x, n){ 
  x <- subset(x, x$V5 > n )
  return(x)
}

############################# read data and process ##############################

# disclaimer- this is the non-optimal way of doing all this processing, but it works so ¯\_(ツ)_/¯

library(data.table) # allows use of fread function 

# reading all this data into memory takes several minutes  
# these are processed unmethylated/methylated 
s1 = fread("results_1_cpg.txt.txt")
s2 = fread("results_2_cpg.txt.txt")
s3 = fread("results_3_cpg.txt.txt")
s4 = fread("results_4_cpg.txt.txt")
s5 = fread("results_5_trunc.txt")
s6 = fread("results_6_trunc.txt")
s7 = fread("results_7_trunc.txt")
s8 = fread("results_8_trunc.txt")

# creates a list of data frames, where the names of each item in the list is the same as the dataframe name 
sample_list <- list(s1 = s1, s2 = s2, s3 = s3, s4 = s4, s5 = s5, s6 = s6, s7 = s7, s8 = s8)

# filter each dataframe, and make a key to merge on (consisting of chr and start concantendated)
# then calculate % methylatedfor each cpg. Then only keep the key and % methylated columns 
sample_list <- lapply(sample_list, function(x) {x <- filter(x, 5)})
sample_list <- lapply(sample_list, function(x) { x <- make.key(x)})
sample_list <- lapply(sample_list, function(x) { x <- percent(x)})
sample_list <- lapply(sample_list, function(x) { x <- keep(x)})

# take the processed list of dataframes and return it to memory 
list2env(sample_list, .GlobalEnv)

################################## plot data ########################################

##### plot a histogram of how many loci there are at each amount of % methylated 

bins <- seq(0, 100) # sets number of bins- one bin for each integer percent

png(filename="percent_meth_histogram.png") # save histogram 

# adjusted color histogram so that the bars are overlayed and slightly transparent (determined by alpha.f parameter)
hist(s1$percent, col = adjustcolor( "red", alpha.f = 0.5), freq=FALSE, breaks=bins, ylim=c(0, 0.08), xlab = "% methylation", ylab="proportion", main="CpG methlyation")
hist(s2$percent, col = adjustcolor( "orange", alpha.f = 0.2), breaks=bins, freq=FALSE, add=T) # overlay on previous 
hist(s3$percent, col = adjustcolor( "yellow", alpha.f = 0.2), breaks=bins, freq=FALSE, add=T)
hist(s4$percent, col = adjustcolor( "green ", alpha.f = 0.2), breaks=bins, freq=FALSE, add=T)
hist(s5$percent, col = adjustcolor( "blue", alpha.f = 0.2), breaks=bins, freq=FALSE, add=T)
hist(s6$percent, col = adjustcolor( "purple", alpha.f = 0.2), breaks=bins, freq=FALSE, add=T)
hist(s7$percent, col = adjustcolor( "pink", alpha.f = 0.2), breaks=bins, freq=FALSE, add=T)
hist(s8$percent, col = adjustcolor( "grey", alpha.f = 0.2), breaks=bins, freq=FALSE, add=T)
legend("topleft", c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"), col=c("Red", "Orange", "Yellow", "Green", "Blue", "Purple", "Pink", "Grey" ), lwd=2)
dev.off()

##### plot a density plot of percent methylation 

# make percent methylated column a dataframe of its own (dfs are easier to plot with ggplot which I find kind of annoying)
vec1 = data.frame(x=s1$percent)
vec2 = data.frame(x=s2$percent)
vec3 = data.frame(x=s3$percent)
vec4 = data.frame(x=s4$percent)
vec5 = data.frame(x=s5$percent)
vec6 = data.frame(x=s6$percent)
vec7 = data.frame(x=s7$percent)
vec8 = data.frame(x=s8$percent)


png("density.png") # save as density.png 
ggplot() + geom_density(data=vec1, aes(x=x, colour="red")) + # plot density for each sample type 
  geom_density(data=vec2, aes(x=x, colour="orange")) + geom_density(data=vec3, aes(x=x, colour="yellow")) + 
  geom_density(data=vec4, aes(x=x, colour="green")) +   geom_density(data=vec5, aes(x=x, colour="blue")) + 
  geom_density(data=vec6, aes(x=x, colour="purple")) +   geom_density(data=vec7, aes(x=x, colour="pink")) + 
  geom_density(data=vec8, aes(x=x, colour="grey")) + 
  scale_color_discrete(name = "sample number", labels = c("ALS1", "ALS2", "ALS3", "ALS4", "CTRL1", "CTRL2", "CTRL3", "CTRL4")) + 
  xlab("% methylation") + theme_light()
dev.off()

################################## merge data ########################################

# the best way I could figure out how to do this is a merge 2 at a time, where the column names are changed each time
# this may take a bit of time 0<time<60 minutes 
merged = merge(s1, s2, by="key", all=TRUE, suffixes = c(".1",".2")) # merging by key 
merged = merge(merged, s3,  by="key", all=TRUE, suffixes = c(".2", ".3"))
merged = merge(merged, s4,  by="key", all=TRUE, suffixes = c(".3", ".4"))
merged = merge(merged, s5,  by="key", all=TRUE, suffixes = c(".4",".5"))
merged = merge(merged, s6,  by="key", all=TRUE, suffixes = c(".5",".6"))
merged = merge(merged, s7,  by="key", all=TRUE, suffixes = c(".6",".7"))
merged = merge(merged, s8,  by="key", all=TRUE, suffixes = c(".7",".8"))

colnames(merged) = c(as.character(seq(0,8))) # change the column names to just the sample number for now 
merged$`0` = NULL # get rid of first column 


################################## pca w/ missing data ########################################

# uses pcamethods library (from bioconductor, see imports)
md <-  prep(merged, scale="none", center=TRUE) 
resPPCA<- pca(md, method="ppca", center=FALSE, nPcs=8) # performs pca, not the best, should be on the covariance 

# plot loadings and pcs 
png(filename="scores-loading.png")
slplot(resPPCA)
dev.off()

png(filename="plotpcs.png")
plotPcs(resPPCA)
dev.off()

# see how good it is predicted from pcs 
predict_resPPCA <- predict(resPPCA, merged)
print(cor(predict_resPPCA$scores[,1], scores(resPPCA))) # see how well the pca performs relative to real data 
 
################################## pca w/o missing data ########################################

# at this point it may be useful to clear the memory and reload just the merged dataset. Up to you. 
# rm(list = ls()) 
# merged = read.table("~/Desktop/zaitlen_lab_desktop/merged_all.txt") 

# remove missing data 
merged <- na.omit(merged) # listwise deletion of missing

# calculate the covariance matrix across the 
cv = cov(na.omit(merged))

# perform PCA 
pca = prcomp(cv, center = TRUE, scale. = TRUE) # uses default R pca method 
plot(pca$x[,1], pca$x[,2]) # plot the first two pcs
summary(pca) # print a summary of the pca, including the % of variance explained by the analysis 


# merged <- scale(merged)
# fit <- kmeans(merged, 8) 
# 
# # chr1 = subset(merged, grepl("chr1", merged$`0`)==TRUE)
# png(filename="scatter.png")
# plot(seq(1, nrow(merged)), merged$`1`, col="blue", pch=18, xlab="genomic position", ylab="% methylation")
# points(seq(1, nrow(merged)), merged$`2`, col="red", pch=18)
# points(seq(1, nrow(merged)), merged$`3`, col="blue", pch=18)
# points(seq(1, nrow(merged)), merged$`4`, col="red", pch=18)
# points(seq(1, nrow(merged)), merged$`5`, col="blue", pch=18)
# points(seq(1, nrow(merged)), merged$`6`, col="red", pch=18)
# points(seq(1, nrow(merged)), merged$`7`, col="blue", pch=18)
# points(seq(1, nrow(merged)), merged$`8`, col="red", pch=18)
# dev.off()

# pca by hand 


