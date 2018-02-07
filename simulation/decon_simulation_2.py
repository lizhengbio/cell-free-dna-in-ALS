# 8 nov 17
# updates to the previous simulation so it makes more sense
# integrated known information about tissues
# author <christa.caggiano@ucsf.edu>

import pandas as pd
import random
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm


def generate_proportion(tissue, individuals):
    """
    generates the contribution of different tissues for
    each individual
    :param tissue: number of tissues to be simulated
    :param individuals: number of individuals
    :return: pandas dataframe of proportion of tissue, of shape individuals by tissue
    """
    rows = []
    for i in range(individuals):
        vals = np.random.choice(10, tissue)  # pick tissue number of values from 1 to 10 to be proportions
        rows.append([x/sum(vals) for x in vals])  # proportions must sum to 1
    return pd.DataFrame.from_records(rows)


def generate_reference(tissue, sites):
    """
    generate reference methylation values for each tissue
    :param tissue: number of tissues
    :param sites: number of CpG sites
    :return: pandas dataframe of proportion of methylated for each site
    """
    rows = []
    for t in range(tissue):  # generate for each tissue
        vals = []
        for s in range(sites):
            vals.append(random.randint(0, 1))
        rows.append(vals)
    return pd.DataFrame(rows).round(decimals=2)

def generate_mixing_proportion(individuals, sites):
    for


if __name__ == "__main__":

    individuals = 5
    sites = 10
    tissues = 2

    # p = generate_proportion(tissues, individuals).as_matrix()
    r = generate_reference(tissues, sites).as_matrix()
    # o = np.around(np.dot(p, r), decimals=2)

    # x, y = o.shape
    # for i in range(x):
    #     for j in range(y):
    #         if random.randrange(0, 1) == 1:
    #             o[i][j] = random.choice([0, 1])
    #         if random.randrange(0, 10) == 1:
    #             o[i][j] = 0.5
    #         if random.randrange(0, 10) == 1:
    #             o[i][j] = 0.75
    #
    # colors = cm.rainbow(np.linspace(0, 1, len(o)))
    # for i, c in zip(range(len(o)), colors):
    #     plt.scatter(range(len(o[i, :])), o[i, :], color=c, label="individual: " + str(i))
    # plt.xlabel("CpG site")
    # plt.ylabel("Proportion methylated")
    # # plt.legend(loc=0)
    # plt.show()
