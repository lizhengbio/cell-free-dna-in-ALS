import matplotlib.pyplot as plt

from optimization.naive_optimization import perform_optimization as naive
from optimization.quadratic_programming_optimization import perform_optimization as qp
from simulation.simulation_utils import *
from utils.error import corr
from utils.run_utils import *

# global simulation parameters
individuals = 1  # number of people optimizing for
sites = 10000  # number of cpg sites
tissues = 5  # number of tissues
read_depth = 100  # read depth (methylated/unmethylated counts)
noise = 0.01


def generate_simulated_optimization(individuals, sites, tissues, read_depth, method, noise):

    proportions = generate_proportion(individuals, tissues).as_matrix()  # initialized proportions of tissue
    print(proportions)
    reference = generate_reference(tissues, sites).as_matrix()  # cpg methylation fraction
    observed = np.matmul(proportions, reference)  # observed estimated is just reference times the proportions

    observed = observed + np.random.normal(0, noise, observed.shape)  # add small amounts of noise to observed
    observed = round_to_one(observed)
    depth = generate_depth(sites, individuals, read_depth)
    methylated = generate_counts(depth, observed, sites, individuals)
    unmethylated = depth - methylated

    proportions_est = np.zeros((individuals, tissues)) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est/(np.sum(proportions_est))

    print(methylated.shape)
    print(unmethylated.shape)
    print(proportions_est.shape)
    print(reference.shape)
    return method(proportions_est, reference, methylated, unmethylated)


def plot_simulation():

    naive_error = []
    qp_error = []

    error = 0
    for individual in range(individuals):
        truth, guess = generate_simulated_optimization(1, sites, tissues, read_depth, naive, noise)
        error += (corr(truth, guess))
    naive_error.append(error/individuals)

    error = 0
    for individual in range(individuals):
        guess = generate_simulated_optimization(1, sites, tissues, read_depth, qp, noise)
        error += (corr(truth, guess))
    qp_error.append(error/individuals)

    # @TODO box and whisker
    plt.plot(noise_range, naive_error, '-bo')
    plt.plot(noise_range, qp_error, '-go')
    plt.ylim(0, 1.2)
    plt.xscale("log")
    plt.show()



