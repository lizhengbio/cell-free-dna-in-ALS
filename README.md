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
    	89* methylation extraction txt files --- used for rest of analyses

### post-analysis of WGBS data
------  

1. Pool any replicates 
2. Sort files 
3. Collapse all observations of a given CpG site 
4. For each CpG site present, calculate percent methylation, i.e. number of methylated observations/total observations

### reference datasets 
------

#### methylation chip datasets

| GEO accession | Study         | File link  | Tissue(s)  |
|---------------|---------------|------------|------------|
| [GSE40279](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE40279)| [Hannum et al (2012)](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/23177740)| [series_matrix](https://goo.gl/CaGDpK)| whole blood|
| [GSE40360](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE40360)   | [Huynh et al (2014)](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/24270187)| [series_matrix](https://goo.gl/vbjtji) | frontal lobe white matter|
| [GSE71378](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE71378) | [Synder et al (2016)](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/26771485)| [series_matrix](https://goo.gl/hvuGWr)     | healthy cfDNA|
| [GSE48472](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE48472) | [Slieker et al (2013)](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-6-26#Bib1)| [series_matrix](https://goo.gl/pP4Zoc)  | various tissues| 
| [GSE62727](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/geo/query/acc.cgi?acc=GSE62727) | [Zhou et al (2017)](https://www-ncbi-nlm-nih-gov.ucsf.idm.oclc.org/pubmed/28849195)| [series_matrix](https://goo.gl/ZbBykZ)  | left atrium| 

#### wgbs and rrbs datasets 

| Study         | File link  |
|---------------|------------|
| [Roadmap Epigenetics Consortium (2015)](www.nature.com/nature/journal/v518/n7539/full/nature14248.html)| [Fractional methylation](http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation.tar.gz), [File key](http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/EG.mnemonics.name.xls) | 






### other 

#### file locations   
parallel job batch submitted to SGE using **qsub_submit.sh**  

original pipeline **/ye/zaitlenlabstore/christacaggiano/dna-me-pipeline/HAIB/bismark_pipeline** 

als files **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/fastq** 

pipeline output **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/FILENAME_HERE**

file key **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/sample_key.txt** 

fastqc scripts **/ye/zaitlenlabstore/christacaggiano/primary_ALS_WGBS/fastqc.txt** 

  
  
