# counts number of lines in final processed bam files
# for ALS data quality control
# fall 2017
# author <christa.caggiano@ucsf.edu>

from __future__ import print_function  # makes it python2 compatible
import glob

# in directory on cluster, get all unsorted bams for counting
files = glob.glob("*_BC_*/unsortedButMerged_ForBismark_file/*_unsorted_merged.bam")

# write number of lines for each file to line_counts.txt
with open("bam_line_counts.txt", "w") as output_file:
    for file in files:
        num_lines = sum(1 for line in open(file))
        print(file+"\t"+str(num_lines), file=output_file)