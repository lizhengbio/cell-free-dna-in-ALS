# hierarchical clustering on top 500 sparse sites 
# jan 2018 
# author <christa.caggiano@ucsf.edu> 

rm(list=ls())
setwd("~/Documents/UCSF_year1/Zaitlen-rotation1/Zaitlen_lab/data/")
sparse = read.table("sparse_sites.txt")

d <- dist(var(sparse), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D2") 
plot(fit, )

