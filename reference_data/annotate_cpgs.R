# adds chromosome location to cpg probes 
# selects tissues of interest in my reference set for ALS experiments 
# jan 2018 
# author <christa.caggiano@ucsf.edu> 

######################### environment ###############################

rm(list=ls())
setwd("~/Documents/UCSF_year1/Zaitlen-rotation1/Zaitlen_lab/reference_data/")

######################### reference data ###############################

illumina = read.table("HumanMethylationSites.txt", header=TRUE, sep=",") # loads in illumina info with probe name and location
obs_cpg = read.table("cpgs.txt", header=TRUE) # a list of cpgs observed in my ALS dataset 
merged_reference = read.table("../data/reference_chip_data_merged.txt", header=TRUE) # my compiled reference cpg probes 

######################### clean data ###############################

illumina_obs = illumina[illumina$cg07881041 %in% obs_cpg$x, ] # remove any probes from illumina data not seen in ALS data
obs_cpg_annot =  subset(merged_reference, merged_reference$X.Sample_geo_accession %in% illumina_obs$cg07881041) # some data not with probe # are in my matrix, so remove those too (mostly points not in a majority of sources) 

 ######################### tissues of interest ###############################

tissues_of_interest = read.table("tissues_of_interest.txt")$V1 # list of geo accesion #s that are experiments on tissues of interest for ALS (muscle, brain, some blood) 
obs_cpg_annot_tissues =  subset(obs_cpg_annot, select=as.vector.factor(tissues_of_interest)) # subset my data to only these tissues 

######################### annotate and save ###############################

obs_cpg_annot_tissues["IlmnID"] = obs_cpg_annot$X.Sample_geo_accession # clean data, make the column names containing cpg probes identical for merging 
illumina_obs["IlmnID"] = illumina_obs$cg07881041 # clean data, make the column names containing cpg probes identical for merging 

cpgs_probe_loc_by_tissue =  merge(obs_cpg_annot_tissues, illumina_obs, by="IlmnID") # merge illumina data with matrix data by probe name
write.table(cpgs_probe_loc_by_tissue, "cpgs_by_tissue.txt") # save 
