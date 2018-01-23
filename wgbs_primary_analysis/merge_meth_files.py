# simple script to concat percent methylation files from indentical samples sequenced in 2 lanes
# dec 2017 
# author <christa.caggiano@ucsf.edu>  

# imports 
import argparse 
import os 

parser = argparse.ArgumentParser()

# takes in fastq files for processing 
parser.add_argument("file_name", type=int, help="file to be processed")
args = parser.parse_args()
file_name = args.file_name

# file that contains file names to be processed 
prefix = open("prefixes4.txt", "r")

# takes the prefix file and turns it in to a more 
file_list = []
for line in prefix: 
    file_list.append(line.rstrip())
name = file_list[file_name] # name of file to be analyzed 

# finds the file using specific file locations
file_1 = (name + "/unsortedButMerged_ForBismark_file/methylation_extraction/" + "CpG_context_" + name + "_unsorted_merged.deduplicated.txt")

# since the 2nd file is a different lane, change lane number 
file_2 = file_1.replace("L001", "L002")

# command to first concatenate files then sort allowing the command to use 80% of the available memory 
cmd = "cat " + file_1 + " " + file_2  + " | sort -S 80%  > " + file_1[0] + "_merged_cpg.txt"

# run command on system 
os.system(cmd)