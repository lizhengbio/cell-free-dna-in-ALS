# calculating percent methylated sites from Bismark output
# author <arya@aryaboudaie.com>
# author <christa.caggiano@ucsf.edu>

from __future__ import print_function, division
from collections import defaultdict
import argparse


def run(input_file, output_file):
    """
    :param input_file: methylation calls -- INPUT MUST BE SORTED
    :param output_file: collapses methylation
    :return: none
    """

    i = 0
    d = defaultdict(lambda: [0,0])

    # for each methylation call, collapse all observations of the particular CpG site
    with open(input_file) as f:
        for line in f:
            if i == 0:
                i += 1
                continue

            # parse the methylation file
            line_split = line.split()
            key = line_split[2]+"|"+line_split[3]

            sign = line_split[1]  # whether the methylation site was present

            if sign == "+":
                d[key][0] += 1
            else:
                d[key][1] += 1
            i += 1

    # output- for each CpG site, chrom location, start, number methylated, number un-methylated
    # percent methylated, and dummy strand info (for homer analysis)
    # @TODO real strand info
    with open(output_file, "w") as out:
        for key, value in d.items():
            chrom, site = key.split("|")
            print("{} {} {} {} {}".format(chrom, site, value[0], value[1], value[0] + value[1],
                                          value[0]/(value[0]+value[1])*100, "+"), file=out)

if __name__ == "__main__":

    # take in arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("file_name", type=int, help="fastq file for processing")
    args = parser.parse_args()
    file_name = args.file_name  # number from SGE of file to be processed
    file_name = file_name - 1  # subtract 1 because python is 0 indexed

    file_list = []
    i = 0
    name = ""

    # open prefixes to find the appropriate file
    with open("prefixes.txt", "r") as f:
        for line in f:
            if i < file_name:
                i += 1
            else:
                name = line.rstrip()
                break

    fn = name  # appropriate file name
    run(fn, "results_{}.txt".format(name))  # call the function

