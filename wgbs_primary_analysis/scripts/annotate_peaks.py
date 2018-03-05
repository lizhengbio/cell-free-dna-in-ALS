# submitting qsub jobs for Homer annotating peaks processing
# uses Homer http://homer.ucsd.edu/homer/index.html
# author <christa.caggiano@ucsf.edu> 
# 5 Oct 2017 

import argparse
import subprocess

# takes arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument("file", type=int)
args = parser.parse_args()

# file to be processed
file = args.file

# open list of fastq files to be annotated
file_list = []
with open("fastq_files.txt", "r") as f:
    for line in f:
        file_list.append(line.rstrip())

# generate homer command
annotate_cmd = "annotatePeaks.pl" + " " + file_list[file] + " " + "hg38 " + "-annStats " + file + \
               "_stats.txt  > " + file_list[file] + "_homer.txt"

# call command, Homer must be installed
subprocess.call(annotate_cmd, shell=True)

