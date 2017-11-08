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
        # half_time = random.randint(0, 1)
        # if half_time == 1:

        vals = np.random.choice(10, tissue)
        rows.append([x/sum(vals) for x in vals])
    return pd.DataFrame.from_records(rows)


def generate_reference(tissue, sites):
    rows = []
    for t in range(tissue):
        vals = []
        for s in range(sites):
            vals.append(random.randint(0,1))
            # half_time = random.randint(0, 1)
            # if half_time == 1:
            #     vals.append(random.randint(0, 1))
            # else:
            #     a = 0.5
            #     b = 0.5
            #     vals.append(np.random.beta(a, b))
        rows.append(vals)
    return pd.DataFrame(rows).round(decimals=2)


if __name__ == "__main__":

    individuals = 5
    sites = 1000
    tissues = 10

    p = generate_proportion(tissues, individuals).as_matrix()
    r = generate_reference(tissues, sites).as_matrix()
    o = np.around(np.dot(p, r), decimals=2)

    x, y = o.shape
    for i in range(x):
        for j in range(y):
            if random.randrange(0, 1) == 1:
                o[i][j] = random.choice([0, 1])
            if random.randrange(0, 10) == 1:
                o[i][j] = 0.5
            if random.randrange(0, 10) == 1:
                o[i][j] = 0.75

    colors = cm.rainbow(np.linspace(0, 1, len(o)))
    for i, c in zip(range(len(o)), colors):
        plt.scatter(range(len(o[i, :])), o[i, :], color=c, label="individual: " + str(i))
    plt.xlabel("CpG site")
    plt.ylabel("Proportion methylated")
    # plt.legend(loc=0)
    plt.show()
