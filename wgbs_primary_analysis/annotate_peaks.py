# submitting qsub jobs for homer annotating peaks processing 
# author <christa.caggiano@ucsf.edu> 
# 5 Oct 2017 

import argparse
import subprocess
import os
import linecache

parser = argparse.ArgumentParser()
parser.add_argument("input", type=int)
args = parser.parse_args()

input = args.input 

f = open("fastq_files.txt", "r")

file_list = []
for line in f: 
	file_list.append(line.rstrip())


annotate_cmd = "annotatePeaks.pl" + " " + file_list[input] + " " + "hg38 " + "-annStats " + input + "_stats.txt  > " file_list[input] + "_homer.txt"
subprocess.call(annotate_cmd, shell=True)

