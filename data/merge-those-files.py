import csv
import os
"""
This program merges two files containing site positions, one with a list of large dmr regions of the genome with many
methylation sites, and one with a list of single methylation sites. The program adds on the averege methylation values
for the tissues to the file with the dmr regions corresponding with those methylation sites.
"""

# New file to be written, delete old one first
os.system("rm file3.csv")

with open("WGBS_DMRs_v2.txt") as dmr_file, open("cpgs_by_tissue.txt") as meth_file, open("file3.csv", "w") as new_csv:
    dmr_regions = csv.reader(dmr_file, delimiter="\t")
    methylation_sites = csv.reader(meth_file, delimiter="\t")
    new_file = csv.writer(new_csv)

    # For ease of remembering which are the small regions and which are the larger ones, I'll refer to them
    # as big and small
    big_line = next(dmr_regions)  # header for dmr file
    small_line = next(methylation_sites)  # header for methlyation sites
    new_file.writerow(big_line + small_line[3:])  # write the combination of the headers

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


    # The strategy is to read through the lines in turn, starting with the first line of each, and advancing them
    # according to a set of rules to find the corresponding methylation sites for each dmr region.
    # Runs in O(n+m) time, where n is number of dmr regions, and m is number of methylation sites.
    finished_through_files = False
    while not finished_through_files:
        try:
            # chromosome of methylation site
            small_chr = chr_int(small_line[0])
            # start and end of methylation site (end = start+1)
            small_start, small_end = flint(small_line[1]), flint(small_line[2])

            # chromosome of dmr
            big_chr = chr_int(big_line[0])
            # smart and end of dmr (could be arbitrarily large)
            big_start, big_end = flint(big_line[1]), flint(big_line[2])

            # If chromosomes don't match, then advance the smaller line.
            if small_chr < big_chr:
                small_line = next(methylation_sites)
                continue
            elif big_chr < small_chr:
                # Every time we advance the DMR file, we must first write it's contents to the new file, since
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
