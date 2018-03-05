#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r y                            
#$ -j y                           
#$ -l mem_free=20G                 
#$ -l arch=linux-x64               
#$ -l netapp=2G,scratch=2G         
#$ -l h_rt=72:00:00      
#$ -t 1

# controls R script to perform hierarchical clustering
# fall 2017
# author <christa.caggiano@ucsf.edu>

R CMD BATCH hierarchical.R