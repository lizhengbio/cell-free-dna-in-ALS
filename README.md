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

| GEO accession | Study         | File link  | Tissue(s)  |
| ------------- |---------------| -----------|------------|

| [GSE40279](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE40279)| [Hannum et al (2012)](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/23177740)| [series_matrix](https://goo.gl/CaGDpK)| whole blood|

| [GSE40360](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE40360)   | [Huynh et al (2014)](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/24270187)| [series_matrix(ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71378/matrix/GSE71378_series_matrix.txt.gz) | frontal lobe white matter|

| [GSE71378](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE71378) | [Synder et al (2016)](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/26771485)| [series_matrix](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40360/matrix/GSE40360_series_matrix.txt.gz)     | healthy cfDNA|

| [GSE48472](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE48472) | [Slieker et al (2013)](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-6-26#Bib1)| [series_matrix](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48472/matrix/GSE48472_series_matrix.txt.gz)  | various tissues| 

### other 

#### file locations   
parallel job batch submitted to SGE using **qsub_submit.sh**  

original pipeline **/ye/zaitlenlabstore/christacaggiano/dna-me-pipeline/HAIB/bismark_pipeline** 

als files **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/fastq** 

pipeline output **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/FILENAME_HERE**

file key **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/sample_key.txt** 

fastqc scripts **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/fastqc.txt** 

  
  
