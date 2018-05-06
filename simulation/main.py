from simulation_utils import *
from error import mse, scaled_mse, corr
from naive_optimization import perform_optimization as naive
from quadratic_programming_optimization import perform_optimization as qp
import numpy as np
import matplotlib.pyplot as plt


def generate_optimization(individuals, sites, tissues, method):

    proportions = generate_proportion(individuals, tissues).as_matrix()  # initialized proportions of tissue
    reference = generate_reference(tissues, sites).as_matrix()  # cpg methylation fraction
    observed = np.matmul(proportions, reference)  # observed estimated is just reference times the proportions

    depth = generate_depth(sites, individuals, read_depth)
    methylated = generate_counts(depth, observed, sites, individuals)
    unmethylated = depth - methylated

    observed = observed + np.random.normal(0, 0.05, observed.shape)  # add small amounts of noise to observed

    proportions_est = np.zeros((individuals, tissues)) + 0.5  # start with random estimation of proportions

    return method(proportions_est, proportions, observed, reference, methylated, unmethylated)

if __name__ == "__main__":

    individuals = 10 # just optimizing for one person
    sites = 10000  # number of cpg sites
    tissues = 100  # number of tissues
    read_depth = 10



    site_range = [10, 100, 1000, 10000, 100000]

    naive_error = []
    qp_error = []

    for site in site_range:
        error = 0
        for individual in range(individuals):
            truth, guess = generate_optimization(1, sites, tissues, naive)
            error += (corr(truth, guess))
        naive_error.append(error/individuals)

    for site in site_range:
        error = 0
        for individual in range(individuals):
            truth, guess = generate_optimization(1, sites, tissues, qp)
            error += (corr(truth, guess))
        qp_error.append(error/individuals)

    print(naive_error)
    print(qp_error)

    plt.plot(site_range, naive_error, '-bo')
    plt.plot(site_range, qp_error, '-go')
    plt.ylim(0, 1)
    plt.xscale("log")
    plt.show()