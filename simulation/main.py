from simulation_utils import *
from error import mse, scaled_mse, corr
from naive_optimization import perform_optimization as naive
from quadratic_programming_optimization import perform_optimization as qp
import numpy as np
import matplotlib.pyplot as plt


def generate_optimization(individuals, sites, tissues, read_depth, method):

    proportions = generate_proportion(individuals, tissues).as_matrix()  # initialized proportions of tissue
    reference = generate_reference(tissues, sites).as_matrix()  # cpg methylation fraction
    observed = np.matmul(proportions, reference)  # observed estimated is just reference times the proportions

    # @TODO check between 0 1
    # observed = observed + np.random.normal(0, 0.05, observed.shape)  # add small amounts of noise to observed

    depth = generate_depth(sites, individuals, read_depth)
    methylated = generate_counts(depth, observed, sites, individuals)
    unmethylated = depth - methylated

    proportions_est = np.zeros((individuals, tissues)) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est/(np.sum(proportions_est))

    return method(proportions_est, proportions, observed, reference, methylated, unmethylated)

if __name__ == "__main__":

    individuals = 100  # number of people optimizing for
    sites = 1000  # number of cpg sites
    tissues = 50  # number of tissues
    read_depth = 100  # read depth (methylated/unmethylated counts)

    tissue_range = [10, 100, 1000]

    naive_error = []
    qp_error = []

    for tissues in tissue_range:
        error = 0
        for individual in range(individuals):
            truth, guess = generate_optimization(1, sites, tissues, read_depth, naive)
            error += (corr(truth, guess))
        naive_error.append(error/individuals)
        print(naive_error)

    for tissues in tissue_range:
        error = 0
        for individual in range(individuals):
            truth, guess = generate_optimization(1, sites, tissues, read_depth, qp)
            error += (corr(truth, guess))
        qp_error.append(error/individuals)
        print(qp_error)

    print(naive_error)
    print(qp_error)

    # box and whisker
    plt.plot(tissue_range, naive_error, '-bo')
    plt.plot(tissue_range, qp_error, '-go')
    plt.ylim(0, 1.2)
    plt.xscale("log")
    plt.show()