# author <christa.caggiano@ucsf.edu>

import pandas as pd
import random
import numpy as np
import scipy
import scipy.stats as stats
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error

def generate_proportion(individuals, tissue):
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
            vals.append(random.uniform(0, 1))
        rows.append(vals)
    return pd.DataFrame(rows)


def log_likelihood(proportions, observed, reference):

    # proportions = proportions/(proportions.sum(axis=0))
    b = np.transpose(np.matmul(proportions, reference))
    sigma = np.var(observed-b)
    N = len(proportions)

    ll = np.sum(-N/2*np.log(2*np.pi) - N/2 * np.log(sigma) - 1/(2*sigma**2) * (np.sum((observed - b)**2)))

    return ll


def log_likelihood2(proportions, observed, reference):

    proportions = proportions/(proportions.sum(axis=0))
    b = np.transpose(np.matmul(proportions, reference))
    sigma = np.var(observed-b)
    N = len(proportions)

    ll = np.sum(-N/2*np.log(2*np.pi) - N/2 * np.log(sigma) - 1/(2*sigma**2) * (np.sum((observed - b)**2)))

    return ll


def perform_optimization_qp(individuals, tissues, proportions_est ):

    bounds = tuple((0, 1) for x in range(np.shape(proportions_est)[1]))
    cons = ({'type': 'eq', 'fun': lambda x: 1 - sum(x)})

    prop_guess_1 = minimize(log_likelihood, proportions_est, args=(observed, reference),
                            bounds=bounds, constraints=cons, method="SLSQP", options={'maxiter': 1000})

    return mean_squared_error(np.transpose(proportions[0]), (prop_guess_1["x"]/np.sum(prop_guess_1["x"])))


def perform_optimization_lame(individuals, tissues, proportions_est):

    bounds = tuple((0, 1) for x in range(np.shape(proportions_est)[1]))

    prop_guess_2 = minimize(log_likelihood2, proportions_est, args=(observed, reference),
                            bounds=bounds, method="L-BFGS-B", options={'maxiter': 1000})

    return mean_squared_error(np.transpose(proportions[0]), (prop_guess_2["x"]/np.sum(prop_guess_2["x"])))



if __name__ == "__main__":

    individuals = 1
    sites = 10
    tissues = 100
    read_depth = 100

    reg = []
    qp = []

    log_scale = [10**i for i in range(1,2)]

    for sites in range(0, 20000, 1000):

        proportions = generate_proportion(individuals, tissues).as_matrix()
        reference = generate_reference(tissues, sites).as_matrix()
        observed = np.matmul(proportions, reference)

        observed = observed + np.random.normal(0, 0.05, observed.shape)

        proportions_est = np.random.rand(individuals, tissues)

        reg.append(perform_optimization_lame(individuals, tissues))
        qp.append(perform_optimization_qp(individuals, tissues))

    print(reg)
    print(qp)


