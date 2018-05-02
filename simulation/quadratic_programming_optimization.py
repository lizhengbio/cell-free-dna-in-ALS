# naive optimization that constrains by summing row and divide by sum
# May 2018
# author <christacaggiano@ucsf.edu>

# imports
import pandas as pd
import random
import numpy as np
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
    """
    calculate log likelihood for optimization

    :param proportions: current estimation of proportions
    :param observed: observed methylation states
    :param reference: reference methylation by tissue
    :return: log likelihood
    """

    b = np.transpose(np.matmul(proportions, reference))
    sigma = np.var(observed-b)
    N = len(proportions)

    # log likelihood function, copied from n. zaitlen's r code
    ll = np.sum(-N/2*np.log(2*np.pi) - N/2 * np.log(sigma) - 1/(2*sigma**2) * (np.sum((observed - b)**2)))

    return ll


def perform_optimization(proportions_est, proportions, observed, reference):
    """
    performs optimization using quadratic programming

    :param proportions_est: estimation of proportions
    :param proportions: 'true' proportions
    :param observed: methylation patterns for cpgs
    :param reference: reference methylation patterns for all tissues
    :return: mean squared error between the true and estimated proportions
    """

    bounds = tuple((0, 1) for x in range(np.shape(proportions_est)[1]))
    cons = ({'type': 'eq', 'fun': lambda x: 1 - sum(x)})

    prop_guess = minimize(log_likelihood, proportions_est, args=(observed, reference),
                            bounds=bounds, constraints=cons, method="SLSQP")

    # return mean_squared_error(np.transpose(proportions[0]), prop_guess["x"])
    return np.transpose(proportions[0]), prop_guess["x"], mean_squared_error(np.transpose(proportions[0]), (prop_guess["x"]))


if __name__ == "__main__":

    individuals = 1  # just optimizing for one person
    sites = 10000  # number of cpg sites
    tissues = 5  # number of tissues
    read_depth = 100

    # site_range = [10, 100, 1000, 100000]

    proportions = generate_proportion(individuals, tissues).as_matrix()  # randomly initialized proportions of tissue for individual
    reference = generate_reference(tissues, sites).as_matrix()  # cpg methylation fraction
    observed = np.matmul(proportions, reference)  # observed estimated is just reference times the proportions

    observed = observed + np.random.normal(0, 0.05, observed.shape)  # add small amounts of noise to observed

    proportions_est = np.random.rand(individuals, tissues)
    proportions_est = proportions_est/proportions_est.sum() # start with random estimation of proportions
    mean_square_error = perform_optimization(proportions_est, proportions, observed, reference)  # perform optimization and return error
    print(mean_square_error[0])
    print(mean_square_error[1])
    print(mean_square_error[2]/(np.sum(np.square(proportions))))
