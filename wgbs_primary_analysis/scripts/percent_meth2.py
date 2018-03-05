# OLD
# calculates percent methylation
# winter 2017
# author <aryab@google.com>
# author <christa.caggiano@ucsf.edu>


from __future__ import print_function, division
from collections import defaultdict
import progressbar

if __name__ == "__main__":

    # calculates percent methylation for a file already processed by percent_meth.py

    i = 0
    d = defaultdict(lambda: [0,0])

    with progressbar.ProgressBar(max_value=129311257+10000) as bar:  # output a progress bar
        with open("results.txt") as f:
            for line in f:
                line_split = line.split()
                key = line_split[0]+"|"+line_split[1]
                d[key][0]+=int(line_split[2])
                d[key][1]+=int(line_split[3])
                i += 1
                if i % 10000 == 0:
                    bar.update(i)

    # output new results
    with open("results_new.txt", "w") as out:
        for key, value in d.items():
            chrom, site = key.split("|")
            print(f"{chrom} {site} {value[0]} {value[1]} {value[0]/(value[0]+value[1])}", file=out)
