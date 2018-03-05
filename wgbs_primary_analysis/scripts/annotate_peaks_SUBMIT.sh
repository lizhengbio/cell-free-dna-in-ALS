#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r y                            
#$ -j y                           
#$ -l mem_free=20G                 
#$ -l arch=linux-x64               
#$ -l netapp=2G,scratch=2G         
#$ -l h_rt=00:29:30                
#$ -t 1:211

# SGE script for controlling peak annotation
# corresponds to python file annotate_peaks.py
# fall 2017
# author <christa.caggiano@ucsf.edu>
                                 
python annotate_peaks.py $SGE_TASK_ID
echo "DONE"


