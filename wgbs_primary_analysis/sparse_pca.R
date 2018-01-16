# script to perform sparse PCA on large methylation dataset 
# based on ReFACTOR method 
# dec 2017 
# author <christa.caggiano@ucsf.edu> <noah.zaitlen@ucsf.edu> 

############################# environment ##################################

# clear working memory 
rm(list=ls())

# set working directory 
# requires "merged_all.txt"
# setwd("~/Desktop/zaitlen_lab_desktop/")

######################## sparePCA function from N. Zaitlen###################### 
sparsePCARef <- function(O, K, t) {

  # O is data matrix 
  # K is the number of components 
  # t is the number of methylation sites to be considered 

  # more information about ReFActor method here: http://www.cs.tau.ac.il/~heran/cozygene/software/refactor.html
  
  num_components = K;
  
  #RUN PCA
  pcs = prcomp(scale(t(O)));
  
  coeff = pcs$rotation
  score = pcs$x
  
  x = score[,1:K]%*%t(coeff[,1:K]);
  An = scale(t(O),center=T,scale=F)
  Bn = scale(x,center=T,scale=F)
  An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
  Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))
  
  
  # Find the distance of each site from its low rank approximation.
  distances = apply((An-Bn)^2,2,sum)^0.5 ;
  dsort = sort(distances,index.return=T);
  ranked_list = dsort$ix
  
  #run PCA on selected sites
  sites = ranked_list[1:t];
  pcs = prcomp(scale(t(O[sites,])));
  score = pcs$x
  
  
  SPCs = t(score[,1:num_components])
  
  return(SPCs)
}


################################### run PCA #######################################
setwd("~/Documents/UCSF_year1/Zaitlen-rotation1/Zaitlen_lab/data/")
# read in txt containing all CpGs found in all 8 samples used in this study 
# this is a large file that may take several minutes to load into R 
# an alternative to read.table() is: 
# require(data.table)
# merged = fread("merged_all.txt")
merged = read.table("merged_all.txt")

# get the variance across samples for each CpG site
# takes pretty long 
getVar = apply(merged[, -1], 1, var)

# sets a threshold for what variance should be kept 
param=1e-4

# gets sites with a variance larger than the threshold
most_var = merged[getVar > param & !is.na(getVar), ]  

# removes nas 
most_var = na.omit(most_var)

# writes table with results for futher analyses 
write.table(most_var, 'merged-most-var.txt')

# runs sparse PCA on dataset, with 8 PCs and 500 sites
pcs = sparsePCARef(most_var, 8, 500)  
write.table(pcs, "pcs.txt") # write table with the output 

pcs = as.matrix(pcs) # just makes sure pcs are in an easy format to play nicely with plotting 

# i calculated variance explained by hand for each pc using formula - standard.deviation(pc)^2 / sum of all pc sdev^2 * 100 

# sets colors 
colors = c("red", "red", "red", "red", "black", "black", "black", "black")

# plots first pc versus 2nd pc as scatterplot 
plot(pcs[1,], pcs[2,], col=colors, xlab="PC1 81% of var", ylab="PC2 15% of var", pch=19)

# adds text with the number of days each person (if relevant) has had ALS 
text(pcs[1,], pcs[2,], labels=c("987", "560", "377", "652", "0", "0", "0", "0"), pos=1)

# 2nd pc versus 3rd 
plot(pcs[1,], pcs[3,], col=colors, xlab="PC2 15% of var", ylab="PC3 2% of var", pch=19)

# adds legend to plot 
legend("topright",legend=c("ALS", "CTRL"), col=c("aquamarine1", "aquamarine4"), lty=1, cex=0.8)


