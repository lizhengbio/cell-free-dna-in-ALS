# finds cpgs that are most different from baseline (other tissues and blood) 
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

############################### environment ########################################## 
# set wd 
rm(list=ls())
setwd("/Users/Christa.Caggiano/Documents/UCSF_year1/Zaitlen-rotation1/ALS_github/data/")

# imports 
require(data.table)

############################### read in data ########################################## 

# read in matrix of all beta values for all tissues in this experiment 
all = read.csv("reference_matrix.txt", header=TRUE, sep="\t")

# read in leukocyte data 
leuk = fread("/Users/Christa.Caggiano/Documents/UCSF_year1/Zaitlen-rotation1/ALS_github/data/plasma_leuk_0.9.txt")  # sites > 90% methylation 
leuk_names = as.list(colnames(leuk))[-1]  # get leukocyte column names 
leuk_names = leuk_names[-9]

# manually set the identifiers for specific tissues 

cereb = c("GSM1283268", "GSM1283269", "GSM1283270", "GSM1283271", "GSM1283272", "GSM1283273")

muscle = c("GSM1179508", "GSM1179509", "GSM1179510", "GSM1179511", "GSM1179512", "GSM1179513")

prefrontal = c("GSM1283235", "GSM1283237", "GSM1283240", "GSM1283242", "GSM1283243")

white_matter = c("GSM992014", "GSM992015", "GSM992016", "GSM992017", "GSM992018", "GSM992019", 
                 "GSM992020", "GSM992021", "GSM992022", "GSM992023", "GSM992024", "GSM992025", 
                 "GSM992026", "GSM992027", "GSM992028", "GSM992029", "GSM992030", "GSM992031", 
                 "GSM992032")

temporal = c("GSM1283245","GSM1283246", "GSM1283247", "GSM1283248", "GSM1283249", "GSM1283251", 
             "GSM1283252", "GSM1283256", "GSM1283257", "GSM1283258")

# all brain tissue 
brain_tissues = read.table("/Users/Christa.Caggiano/Documents/UCSF_year1/Zaitlen-rotation1/ALS_github/data/braingeo.txt")
brain_tiss = c(levels(brain_tissues$V1))

# manually set 'other tissues' for brain (i.e. tissue that isn't brain)
other_for_brain = read.table("/Users/Christa.Caggiano/Documents/UCSF_year1/Zaitlen-rotation1/ALS_github/data/other_brain.txt")
other_for_brain = c(levels(other_for_brain$V2))


############################### subset by tissue ########################################## 

# make dataframes of only the identifiers relevant to each tissue 
white = all[, white_matter]
brain_tiss = all[, brain_tiss]
cereb  = all[, cereb]
musc = all[, muscle]
pref = all[, prefrontal]
temporal = all[, temporal]
other = all[, other_for_brain]
leukocytes = all[, unlist(leuk_names)]

# calculate the average methylation value for all replicates of a given tissue 
all$leuk_average = rowMeans(leukocytes, na.rm=TRUE)
all$all_brain_average = rowMeans(brain_tiss, na.rm=TRUE)
all$cereb_average = rowMeans(cereb, na.rm=TRUE)
all$musc_average = rowMeans(musc, na.rm=TRUE)
all$pref_average = rowMeans(pref, na.rm=TRUE)
all$other_brain_avg = rowMeans(other, na.rm=TRUE)
all$white_matter = rowMeans(white, na.rm=TRUE)
all$temporal = rowMeans(temporal, na.rm=TRUE)


###############################  top brain  ########################################## 

# subset to relevant brain tissues/other brain tissues and metadata columns 
averages = c("cereb_average", "all_brain_average", "other_brain_avg", "pref_average", "leuk_average", "white_matter", "temporal", "chr", "start", "X.Sample_geo_accession")
brain = all[, averages]

brain = subset(brain, brain$leuk_average > 0.9) # all leukocytes greater than 90% 
brain = subset(brain, brain$all_brain_average < 0.5) # white matter cpgs that are less than 50% methylated 
brain = subset(brain, brain$white_matter > 0.5) # white matter cpgs that are less than 50% methylated 
brain = subset(brain, brain$other_brain_avg > 0.8) # other brain tissues must be greater than 80% 

write.table(brain, "brain_cpgs.csv", sep=",", row.names = F)  # write output as csv 

###############################  top muscle  ########################################## 

# other tissues for muscle 
other_columns = c("chr", "start", "end", "leuk_average", "X.Sample_geo_accession", 
                  "X", "gene_id", "closest.SNP", "SNP.allele.freq", "gene_name")

# tissues that are not muslce or leukocytes for comparison 
not_musc_not_leuk = subset(colnames(all), !colnames(all) %in% c(leuk_names, other_columns) & !colnames(all) %in% musc)
not_musc = all[, not_musc_not_leuk]

# average of not muscle tissues 
all$not_musc_avg= rowMeans(not_musc)

# subset to relevant tissues for muscle analysis and metadata 
averages = c("not_musc_avg", "musc_average", "leuk_average", "chr", "start", 
             "X.Sample_geo_accession", "gene_id", "closest.SNP", "SNP.allele.freq", "gene_name")

# subset to muscles 
musc_cpgs = all[, averages]
musc_cpgs = subset(musc_cpgs, musc_cpgs$leuk_average > 0.9)  # all leukocytes over 90% 
musc_cpgs = subset(musc_cpgs, musc_cpgs$musc_average < 0.5) # all muscle cpgs less than 50% methylated 
musc_cpgs = subset(musc_cpgs, musc_cpgs$not_musc_avg > 0.8) # all other tissues must be more than 80% methylated 

write.table(musc_cpgs, "musc_cpgs.csv", sep=",", row.names = F) # write output 

############################### leukocytes  ########################################## 

not_leuk = subset(colnames(all), !colnames(all) %in% c(leuk_names, other_columns)) # all tissues that are not leukocytes 
all$not_leuk_avg = rowMeans(all[, not_leuk]) # all tissues other than leukocytes 
all$leuk_dist = abs(all$leuk_average-all$not_leuk_avg)  # distance between average leukocytes and othe

averages = c("not_leuk_avg", "leuk_dist",
             "leuk_average", "chr", "start", "X.Sample_geo_accession",
             "gene_id", "closest.SNP", "SNP.allele.freq", "gene_name")
leuk_avg = all[, averages]
leuk_avg = all[, averages]

leuk_avg = leuk_avg[rev(order(leuk_avg$leuk_dist)),]
leuk_avg = leuk_avg[complete.cases(leuk_avg), ]
write.table(leuk_avg, "leuk_dist.csv", sep=",", row.names = F)
write.table(a, "leuk_all", sep=",", row.names = F)

# leuk_avg = subset(leuk_avg, leuk_avg$leuk_average <0.5)
# leuk_avg = subset(leuk_avg, leuk_avg$not_leuk_avg > 0.9)
