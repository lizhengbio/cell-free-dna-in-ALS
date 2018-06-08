# gathering some statistics/information about methylation data arrays 
# for python simulation 
# 30 aug 17
# author <christa.caggiano@ucsf.edu> 

rm(list = ls())

# data mentioned in refactor paper, accession number GSE35069, downloaded http manually 
real_data = read.table("/Users/christacaggiano/Documents/UCSF year 1/Zaitlen-rotation1/GSE35069_Matrix_signal_intensities.txt", sep="\t", header=TRUE)

# selects only methylation info 
methylation_data = real_data[ , c(FALSE,FALSE, TRUE) ]

# takes the mean of all samples 
column_means = colMeans(methylation_data)

# reports statistics about each mean 
mean(column_means)
median(column_means)
sd(column_means)

max_value = max(apply(methylation_data, 2, function(x) max(x, na.rm = TRUE))) 
