# given a bed file with cpg count observations,
# pull out cpgs of interest 
# march 2018 
# author <christa.caggiano@ucsf.edu> 

# imports
import csv
import argparse
from pybedtools import BedTool  # requires pybedtools to be installed

def add_range(search_range, file_name):
    """
    for an input bed file containing cpgs of interest, find a user-given range around that cpg to
    search for
    :param search_range: int of search range around the cpg of interest
    :param file_name: bed file name
    :return new_file_name: name of new bed file that represents a range around each cpgs of interest
    """

    new_file_name = make_new_file_name(file_name, search_range)  # generate new file name

    # open the bed file of interest and the output file to be generated
    with open(file_name) as f, open(new_file_name, "w") as out:

        out_bed = csv.writer(out, delimiter="\t")  # bed files must be tab delimited 

        for line in f:

            line = line.split()

            if line[0] == "chr":
                out_bed.writerow(line)  # write the header

            else:

                new_start = int(line[1]) - search_range  # make new search range
                new_end = int(line[2]) + search_range

                new_line = [line[0], new_start, new_end]+line[3:]

                out_bed.writerow(new_line)  # write new line to output file

    return new_file_name


def make_new_file_name(file_name, search_range):
    """
    makes the new file name for output given a few scenarios of input names
    :param file_name: input file name
    :param search_range: the number of bps to inspect around a cpg
    :return: new file name
    """

    # makes new file name if the original file contained .bed
    if ".bed" in file_name:
        file_name = file_name.replace(".bed", "")

    # make new file name if the original file is just text
    if ".txt" in file_name:
        file_name = file_name.replace(".txt", "")

    # output file name is a bed file
    return file_name + "_" + str(search_range) + ".bed"


def bed_intersect(cpg, data):
    """
    using pybedtools, perform an intersection on our cpg file of interest and cpg count observation file
    :param cpg: cpg file
    :param data: count file
    :return: None
    """
    output_name = cpg.replace(".bed", "_results.txt")  # make results file

    # make bedtool objects out of cpg and data files
    a = BedTool(cpg)
    b = BedTool(data)

    # perform a left outer join intersect and move to an output file
    a.intersect(b, loj=True).moveto(output_name)


if __name__ == "__main__":

    # take command line arguments
    parser = argparse.ArgumentParser(description="takes in a bed file with coordinates to search ALS cfDNA data")
    parser.add_argument("cpg_file_name", help="bed file of interesting cpgs", type=str)
    parser.add_argument("data_name", help="bed file of data to search", type=str)
    parser.add_argument("range", help="+/- values to search around bed coordinates", type=int)

    args = parser.parse_args()

    cpg = args.cpg_file_name  # bed file of cpgs of interest
    search_range = args.range  # amount of bp to search around given cpg
    data = args.data_name  # observational cpg counts for als data

    new_file_name = add_range(search_range, cpg)  # make new bed file with a range

    bed_intersect(new_file_name, data)  # perform intersect


