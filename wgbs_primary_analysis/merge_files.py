# merge files 

import argparse 
import os 

parser = argparse.ArgumentParser()

parser.add_argument("file_name", type=int, help="fastq file for processing")
args = parser.parse_args()
file_name = args.file_name
index = file_name - 1

prefix = open("prefixes4.txt", "r")

file_list = []
for line in prefix: 
    file_list.append(line.rstrip())
name = file_list[index]

file_1 = (name + "/unsortedButMerged_ForBismark_file/methylation_extraction/" + "CpG_context_" + name + "_unsorted_merged.deduplicated.txt")
file_2 = file_1.replace("L001", "L002")

cmd = "cat " + file_1 + " " + file_2  + " | sort -S 80%  > " + file_1[0] + "_merged_cpg.txt"

os.system(cmd)
