# 29 aug 17
# simple script to produce simulated observed and reference data
# for given mixtures of origin tissue
# author <christa.caggiano@ucsf.edu>

import numpy as np
import argparse
import random


def choose_random_cases(case):

    """
    pseudo-randomly decide if a case should be selected
    :param case: number of cases (individual or cpgs)
    :return: a dict that contains a bool indicating whether or not case is selected
    """
    case_dict = {i: random.choice([False, True]) for i in range(case)}
    return case_dict


def create_array(case_dict, mean, sd, max, rows, columns, distribution_type, scaling, column_stack):
    """
    creates an array where if a particular case is selected, values are generated around a given distribution
    otherwise values are generated so that they are similar, with a small degree of noise
    :param case_dict: either cpg or individual dictionaries that indicate which cases should be selected
    :param mean:
    :param sd:
    :param max:
    :param size: number of values to be generated
    :param distribution_type:
    :param column_stack: bool value indicating how to appropriately arrange the data
    :return: an array containing values (methylation or cell type proportions) that have a particular distribution
    """
    if column_stack:
        array_list = [[] for i in range(columns)]
        for case in range(len(array_list)):
            if case_dict[case]:
                array_list[case] = distribution_type(mean, sd, rows) * scaling
            else:
                static_value = max / 100 * random.randint(10, 99)
                array_list[case] = [
                    static_value + random.choice(np.arange(-(static_value / 10), static_value / 10, scaling / 10))
                    for i in range(rows)]
            print(array_list)

        return np.column_stack(array_list)

    else:
        array_list = [[] for i in range(rows)]

        for case in range(len(array_list)):
            if case_dict[case] == 1:
                array_list[case] = distribution_type(mean, sd, columns)*scaling
            else:
                static_value = max / 100 * random.randint(10, 99)
                array_list[case] = [static_value + random.choice(np.arange(-(static_value/10), static_value/10, scaling/10))
                                    for i in range(columns)]
        return np.vstack(array_list)

if __name__ == "__main__":

    # accepts user input
    # all choices must be greater than 1
    # default is 5 individuals, 3 tissues, and 3 CpGs
    # parser = argparse.ArgumentParser()
    # parser.add_argument("individuals", help="number of individuals", type=int, choices=range(2,100), default=5)
    # parser.add_argument("tissues", help="number of tissues", type=int, choices=range(2,100), default=3)
    # parser.add_argument("cpgs", help="number of CpG sites", type=int, choices=range(2,100), default=3)

    # args = parser.parse_args()

    # individual = args.individuals
    # cpg = args.cpgs
    # tissue = args.tissues

    individual = 3
    cpg = 4
    tissue = 3
    dmr_num = 0
    control_num = 0


    # TODO proportion methylated
    mean = 7000
    maximum = 60000
    sd = 435

    cpg_dict = choose_random_cases(cpg)  # select random cpgs to be differentially methylated
    individual_dict = choose_random_cases(individual)  # select individuals to have non-normal levels of cfDNA


    # TODO randomly pick one tissue to be diff methylated
    # creates a reference array where methylation values are normally distributed if they are differentially methylated
    # and base these methylation values around the descriptive statistics of the aforementioned data set
    ref_array = create_array(cpg_dict, mean, sd, maximum, tissue, cpg, np.random.normal, 1, True)

    # create an array of cell proportions for each individual where proportions have a lognormal(skewed) distribution
    # indicating a non-normal (perhaps disease) state
    cell_proportion_array = create_array(individual_dict, 2, 0.2, 1, individual, tissue, np.random.lognormal, 0.1, False)
    print(ref_array)
    print(cell_proportion_array)
    observed_array = np.dot(cell_proportion_array, ref_array)
    # print(observed_array)
    # TODO add in noise y = mr + e

