# calculating percent methylated sites
# author <christa.caggiano@ucsf.edu>
# author <aryab@google.com>
from __future__ import print_function, division
from collections import defaultdict
import glob

 
# glob.glob("SRR*/unsortedButMerged_ForBismark_file/methylation_extraction/CpG_context_SRR*_unsorted_merged.deduplicated.txt")
def write_results(key, value, output_file):
    with open(output_file, "a") as out:
        print("{} {} {} {} {}".format(key[0], key[1], value[0], value[1], value[0]/(value[0]+value[1])), file=out)


def run(input_file, output_file):
    i = 0
    previous_key = []
    pluses = 0
    minuses = 0
    with open(input_file) as f:
        for line in f:
            if i == 0:
                i += 1
                continue
            line_split = line.split()
            key = [line_split[2], line_split[3]]
            if previous_key and key != previous_key:
                write_results(previous_key, [pluses, minuses], output_file)
                pluses = 0
                minuses = 0

            previous_key = key

            sign = line_split[1]
            if sign == "+":
                pluses += 1
            else:
                minuses += 1
            i += 1


if __name__ == "__main__":
    filenames = glob.glob("SRR*/unsortedButMerged_ForBismark_file/methylation_extraction/CpG_context_SRR*_unsorted_merged.deduplicated.txt")
    for fn in filenames:
        run(fn, fn.split("/CpG_context_SRR")[0]+"/results.txt")
