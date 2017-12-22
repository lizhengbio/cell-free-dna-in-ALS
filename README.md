# Guide to BS-seq Bioinformatics Workflow
## Zaitlen ALS Project Fall 2017 
## author: <christa.caggiano@ucsf.edu> 

### Whole Genome Bisulfite Sequencing (WGBS) Analysis 
------
Pipeline to process raw WGBS here is based on the WGBS pipeline from the ENCODE consortium, [found on github](https://github.com/ENCODE-DCC/dna-me-pipeline/blob/master/HAIB/bismark_pipeline/bismark_pipeline_main.sh) 

My implementation of this pipeline can be found [here](https://github.com/christacaggiano/ENCODE_WGBS) and uses the following steps: 

0. Check quality of raw fastq files with fastqc 
1. Remove adapters and trim reads
2. Align to hg38 using Bowtie2  
3. Find methylated and unmethylated sites using Bismark methylation extraction  

Output files: 
    * aligned bams 
    * methylation reports (command, summary of run)  
    * methylation extraction txt files --- used for rest of analyses

### post-analysis of WGBS data
------  

1. Pool any replicates 

2. Sort files 
3. Collapse all observations of a given CpG site 
4. For each CpG site present, calculate percent methylation, i.e. number of methylated observations/total observations

### methylation chip datasets 
------

[GSE40279](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE40279) [Hannum_et_all](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/23177740)


| GEO accession | Study         | File link  |
| ------------- |:-------------:| ----------:|
| [GSE40279](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE40279)      | [Hannum_et_al](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/23177740)| [series_matrix](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz)      |
| col 2 is      | centered      |   $12      |
| zebra stripes | are neat      |    $1      |


### other 

#### file locations   
parallel job batch submitted to SGE using **qsub_submit.sh**  

original pipeline **/ye/zaitlenlabstore/christacaggiano/dna-me-pipeline/HAIB/bismark_pipeline** 

als files **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/fastq** 

pipeline output **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/FILENAME_HERE**

file key **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/sample_key.txt** 

fastqc scripts **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/fastqc.txt** 

  
  
