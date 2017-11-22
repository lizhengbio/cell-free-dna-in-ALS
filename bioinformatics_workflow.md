# Guide to BS-seq Bioinformatics Workflow
## Zaitlen ALS Project Fall 2017 
## author: <christa.caggiano@ucsf.edu> 

### pipeline 
based on https://github.com/ENCODE-DCC/dna-me-pipeline/blob/master/HAIB/bismark_pipeline/bismark_pipeline_main.sh
currently called wgbs_last_step.py - will be updated TBA 
1. Trim-galore: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
2. Run bismark alignment on paired-end trimmed fastqs using Bowtie 2 
3. Bismark CpG methylation extraction 
4. Take Bismark methylation files:
    * sort using bash 
    * pool all methylation calls for a given locus- custom python script 
    * pool by sample if necessary - custom python script 
5. output files: 
    * aligned bams 
    * methylation reports (command, summary of run) 
    * methylation extraction txt files --- used for rest of analyses 

### misc files  
parallel job batch submitted to SGE using **qsub_submit.sh**  

original pipeline **/ye/zaitlenlabstore/christacaggiano/dna-me-pipeline/HAIB/bismark_pipeline** 

als files **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/fastq** 

pipeline output **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/FILENAME_HERE**

file key **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/sample_key.txt** 

  
  
