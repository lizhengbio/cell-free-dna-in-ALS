import csv
import os
"""
This program merges two files containing site positions, one with a list of large dmr regions of the genome with many
methylation sites, and one with a list of single methylation sites. The program adds on the averege methylation values
for the tissues to the file with the dmr regions corresponding with those methylation sites.
"""
dmr_regions = csv.reader(open("WGBS_DMRs_v2.txt"), delimiter="\t")
methylation_sites = csv.reader(open("cpgs_by_tissue.txt"), delimiter="\t")

# New file to be written
os.system("rm file3.csv")
new_file = csv.writer(open("file3.csv", "w"))

big_line = next(dmr_regions)  # header
small_line = next(methylation_sites)  # header
new_file.writerow(big_line + small_line[3:])  # write the header

big_line = next(dmr_regions)
small_line = next(methylation_sites)

prev_row = []


def flint(x):
    """
    Returns float(int(x)) - since the spreadsheet had numbers as floats
    :param x: String representing number
    :return: int of that number
    """
    return int(float(x))


def chr_int(x):
    """
    Turns chromosome number into number, needed since chromosome could be x or y
    :param x: 
    :return: 
    """
    try:
        return int(x)
    except ValueError:
        if x.lower() == "y":
            return 101
        elif x.lower() == "x":
            return 100
        else:
            raise ValueError("Not sure what to do with value {}".format(x))


finished_through_files = False
while not finished_through_files:
    try:
        small_chr = chr_int(small_line[0])
        small_start, small_end = flint(small_line[1]), flint(small_line[2])
        big_chr = chr_int(big_line[0])
        big_start, big_end = flint(big_line[1]), flint(big_line[2])
        # import pdb;pdb.set_trace()
        if small_chr < big_chr:
            small_line = next(methylation_sites)
            continue
        elif big_chr < small_chr:
            new_file.writerow(big_line)
            prev_row = big_line
            big_line = next(dmr_regions)
            continue
        # chr is equal on both

        if small_start < big_start:
            small_line = next(methylation_sites)
            continue
        elif big_start <= small_start <= big_end and big_start <= small_end <= big_end:
            row = big_line + small_line[3:]
            keep_going = True
            while keep_going:
                small_line = next(methylation_sites)
                small_chr = chr_int(small_line[0])
                small_start, small_end = flint(small_line[1]), flint(small_line[2])

                if small_chr != big_chr:
                    keep_going = False
                elif big_start <= small_start <= big_end and big_start <= small_end <= big_end:
                    row += small_line[3:]
                else:
                    keep_going = False

            new_file.writerow(row)
            prev_row = row
            big_line = next(dmr_regions)
            # small line already advanced in the above while loop
            continue
        else:  # small interval not in range, but next small number could be in range
            new_file.writerow(big_line)
            prev_row = big_line
            big_line = next(dmr_regions)
            continue
    except StopIteration:  # Either no more small intervals or no more big intervals
        if prev_row != big_line:
            new_file.writerow(big_line)
        for line in dmr_regions:  # finish writing the big_interval lines
            new_file.writerow(line)
        finished_through_files = True
