from simulation_utils import *
from error import mse, scaled_mse, corr
from naive_optimization import perform_optimization as naive
from quadratic_programming_optimization import perform_optimization as qp
import numpy as np
import matplotlib.pyplot as plt


def generate_optimization(individuals, sites, tissues, read_depth, method, noise):

    # proportions = generate_proportion(individuals, tissues).as_matrix()  # initialized proportions of tissue
    # reference = generate_reference(tissues, sites).as_matrix()  # cpg methylation fraction
    # observed = np.matmul(proportions, reference)  # observed estimated is just reference times the proportions
    #
    # # @TODO check between 0 1
    # observed = observed + np.random.normal(0, noise, observed.shape)  # add small amounts of noise to observed
    # observed = round_to_one(observed)
    # depth = generate_depth(sites, individuals, read_depth)
    # methylated = generate_counts(depth, observed, sites, individuals)
    # unmethylated = depth - methylated

    methylated = np.genfromtxt('obs1_meth.csv', delimiter=',')
    unmethylated = np.genfromtxt('obs1_unmeth.csv', delimiter=',')

    reference = np.genfromtxt('obs1_reference.csv', delimiter=',')

    proportions_est = np.zeros((individuals, 119)) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est/(np.sum(proportions_est))

    return method(proportions_est, reference, methylated, unmethylated)

if __name__ == "__main__":

    individuals = 50  # number of people optimizing for
    sites = 10000  # number of cpg sites
    tissues = 50  # number of tissues
    read_depth = 100  # read depth (methylated/unmethylated counts)

    print(generate_optimization(1, sites, tissues, read_depth, naive, 1))

    # read_depth_range = [1000]
    # noise_range = [0.01, 0.05, 0.1, 0.5]
    #
    # naive_error = []
    # qp_error = []
    #
    # for noise in noise_range:
    #     error = 0
    #     for individual in range(individuals):
    #         truth, guess = generate_optimization(1, sites, tissues, read_depth, naive, noise)
    #         error += (corr(truth, guess))
    #     naive_error.append(error/individuals)
    #     print(naive_error)
    #
    # for noise in noise_range:
    #     error = 0
    #     for individual in range(individuals):
    #         truth, guess = generate_optimization(1, sites, tissues, read_depth, qp, noise)
    #         error += (corr(truth, guess))
    #     qp_error.append(error/individuals)
    #     print(qp_error)
    #
    # print(naive_error)
    # print(qp_error)
    #
    # # box and whisker
    # plt.plot(noise_range, naive_error, '-bo')
    # plt.plot(noise_range, qp_error, '-go')
    # plt.ylim(0, 1.2)
    # plt.xscale("log")
    # plt.show()