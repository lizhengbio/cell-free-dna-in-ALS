from __future__ import print_function
import glob

files = glob.glob("*_BC_*/unsortedButMerged_ForBismark_file/*_unsorted_merged.bam")


with open("linecounts.txt", "w") as output_file:
	for file in files:
		num_lines = sum(1 for line in open(file))
		print(file+"\t"+str(num_lines), file=output_file)