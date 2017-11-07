# 1 nov 17
# updates to the previous simulation so it makes more sense
# author <christa.caggiano@ucsf.edu>

import pandas as pd
import random
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm


def generate_proportion(tissue, individuals):
    rows = []
    for i in range(individuals):
        vals = np.random.choice(100, tissue)
        rows.append([x/sum(vals) for x in vals])
    return pd.DataFrame.from_records(rows)


def generate_reference(tissue, sites, mu_max, sigma):
    a, b = random.sample(range(1, 5), 2)
    rows = []
    for t in range(tissue):
        vals = np.random.beta(a, b, sites)
        rows.append([x for x in vals])
    return pd.DataFrame.from_records(rows)


if __name__ == "__main__":

    individuals = 5
    sites = 10
    tissues = 10

    mu_max = 0.5
    sigma = 1

    p = generate_proportion(tissues, individuals).as_matrix()
    r = generate_reference(tissues, sites, mu_max, sigma).as_matrix()
    o = np.dot(p, r)

    # print(o)
    # colors = cm.rainbow(np.linspace(0, 1, len(o)))
    # for i, c in zip(range(len(o)), colors):
    #     plt.scatter(range(len(o[i, :])), o[i, :], color=c, label="individual: " + str(i))
    # plt.xlabel("CpG site")
    # plt.ylabel("Proportion methylated")
    # # plt.legend(loc=0)
    # plt.show()
