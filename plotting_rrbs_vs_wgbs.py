# plotting information about RRBS and WGBS
# 10 Oct 2017
# author <christa.caggiano@ucsf.edu>

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import numpy as np
from scipy import stats

if __name__ == "__main__":

    sites = pd.read_table("/Users/Christa.Caggiano/Desktop/depth.txt")

    feat = []
    for row in (sites.ix[:, 0]):
        for word in row.split():
            if 'wg' not in word and "total" not in word:
                feat.append(int(word))

    feat_array = np.asanyarray(feat[0:len(feat)-1])
    print("std: {}".format(np.std(feat_array)))
    print("mean: {} ".format((np.mean(feat_array))))
    print("minimum: {}".format(np.min(feat_array)))
    print("maximum: {}".format(np.max(feat_array)))

#     files = glob.glob("*_stats.txt")
#     f = files[0]
#
#     test = pd.read_table(f)
#     x = list(test.loc[:11, "Annotation"])
#     y = test.loc[:11, "Number of peaks"].astype(float)
#
#     for file in files[1:]:
#         test = pd.read_table(file)
#         y = y+test.loc[:11, "Number of peaks"].astype(float)
#
#     y = y/len(files)
#     sns.barplot(x, y, palette="GnBu_d")
#
# plt.savefig('foo.png')
