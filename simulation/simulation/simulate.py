import matplotlib.pyplot as plt

from optimization.naive_optimization import perform_optimization as naive
from optimization.quadratic_programming_optimization import perform_optimization as qp
from simulation.simulation_utils import *
from utils.error import corr, mse, scaled_mse
from utils.run_utils import *


def generate_simulated_optimization(individuals, sites, tissues, read_depth, method, noise, fixed_proportion=None):

    proportions = np.asarray(generate_proportion_fixed(individuals, tissues, fixed_proportion))
    # proportions = generate_proportion(individuals, tissues).as_matrix()  # initialized proportions of tissue
    # print((proportions))
    reference = generate_reference(tissues, sites).as_matrix()  # cpg methylation fraction
    observed = np.matmul(proportions, reference)  # observed estimated is just reference times the proportions

    observed = observed + np.random.normal(0, noise, observed.shape)  # add small amounts of noise to observed
    observed = round_to_one(observed)
    depth = generate_depth(sites, individuals, read_depth)
    methylated = generate_counts(depth, observed, sites, individuals)
    unmethylated = depth - methylated

    proportions_est = np.zeros((individuals, tissues)) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est/(np.sum(proportions_est))

    return proportions, method(proportions_est, reference, methylated, unmethylated)


# global simulation parameters
individuals = 2  # number of people optimizing for
sites = 10000  # number of cpg sites
tissues = 5  # number of tissues
read_depth = 100  # read depth (methylated/unmethylated counts)
noise = 0.01

naive_error = []
qp_error = []

read_depth_range = [10, 100, 1000, 10000]
error = 0
for read_depth in [10, 100, 1000, 10000]:
    for individual in range(individuals):
        truth, guess = generate_simulated_optimization(1, sites, tissues, read_depth, qp, noise, 100)
        error += (scaled_mse(truth[0], guess))
    qp_error.append(error / individuals)
print(qp_error)

qp_error = []

for read_depth in [10, 100, 1000, 10000]:
    for individual in range(individuals):
        truth, guess = generate_simulated_optimization(1, sites, tissues, read_depth, qp, noise, 10)
        error += (corr(truth[0], guess))
    qp_error.append(error / individuals)
print(qp_error)

qp_error = []

for read_depth in [10, 100, 1000, 10000]:
    for individual in range(individuals):
        truth, guess = generate_simulated_optimization(1, sites, tissues, read_depth, qp, noise, 1)
        error += (corr(truth, guess))
    qp_error.append(error / individuals)
print(qp_error)

qp_error = []

for read_depth in [10, 100, 1000, 10000]:
    for individual in range(individuals):
        truth, guess = generate_simulated_optimization(1, sites, tissues, read_depth, qp, noise, .01)
        error += (corr(truth, guess))
    qp_error.append(error / individuals)
print(qp_error)

# print(qp_error)
# # # @TODO box and whisker
# # plt.plot(noise_range, naive_error, '-bo')
# plt.plot([10, 100, 1000, 10000], qp_error, '-go')
# plt.ylim(0, 10)
# plt.xscale("log")
# plt.show()