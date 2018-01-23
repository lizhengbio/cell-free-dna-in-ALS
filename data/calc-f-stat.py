import numpy as np
import pandas as pd
import scipy.stats as stats


def get_f_oneway(als, ctrl):
    return stats.f_oneway(als, ctrl)


def get_rank_sum(als, ctrl):
    return stats.ranksums(als, ctrl)

if __name__=="__main__":

    df = pd.read_table("merged-sites-all.txt", delimiter=" ")

    f_stat = []
    f_p_val = []
    rank_sum = []
    rank_p_val = []

    for index, row in df.iterrows():
        als = np.array(row[1:5])
        ctrl = np.array(row[5:])

        f, p = get_f_oneway(als, ctrl)

        f_stat.append(f)
        f_p_val.append(p)

        r, p = get_rank_sum(als, ctrl)

        rank_sum.append(r)
        rank_p_val.append(p)


    df["fstat"] = f_stat
    df["f-pval"] = f_p_val
    df["rank"] = rank_sum
    df["rank_p_val"] = rank_p_val

    df.to_csv("fstat-merged-sites.csv")