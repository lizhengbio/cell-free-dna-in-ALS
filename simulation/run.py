from optimization.quadratic_programming_optimization import perform_optimization as qp
from utils.run_utils import *


def generate_optimization(reference, methylated, unmethylated, method):

    proportions_est = np.zeros((1, reference.shape[0])) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est / (np.sum(proportions_est))
    return method(proportions_est, reference, methylated, unmethylated)


if __name__ == "__main__":


    for patient in ["ctrl1_common_header.txt", "als1_common_header.txt"]:
        reference, methylated, unmethylated = generate_matrices("/Users/Christa.Caggiano/Desktop/zaitlen_lab_desktop/" + patient, "tissues.txt")
        x = generate_optimization(reference, methylated, unmethylated, qp)
        np.savetxt(patient + "_qp__no_WGBS_results.txt", x)