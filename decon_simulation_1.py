# 29 aug 17
# simple script to produce simulated observed and reference data
# for given mixtures of origin tissue
# author <christa.caggiano@ucsf.edu>

import numpy as np
import argparse
import random

def choose_dmr(cpg):

    """
    pseudo-randomly decide if the cpg should be differentially methylated between tissues
    :param cpg: number of cpgs
    :return: a dict that contains a bool (0 or 1) indicating whether that cpg should be differentially methylated
    """

    cpg_dict = dict.fromkeys(range(cpg))  # generate dict where the keys are the cpgs of interest

    for cpg in cpg_dict:  # for each cpg pick a 0 or 1 and assign it as the value
        cpg_dict[cpg] = random.randint(0, 1)

    return(cpg_dict)


def assign_cpg_to_tissue(cpg_dict, tissue):
    """
    pseudo-randomly assigns a dm cpg to be associated with a particular tissue
    :param cpg_dict: dictionary containing information on whether cpg is differentially methylated
    :param tissue: number of tissues for simulation
    :return: a dict containing which cpgs are associated with a given tissue (for reference simulation)
    """

    # generate tissue data structs
    tissues = list(range(tissue))  # list of each tissue
    tissue_dict = dict.fromkeys(tissues) # makes a blank list for each tissue
    for tissue in tissue_dict:
        tissue_dict[tissue] = []

    # for each cpg if it is a dmr, pick a random tissue and assign that cpg to that tissue key in the tissue dict
    # multiple cpgs may be assigned to one tissue
    for cpg in cpg_dict:
        if cpg_dict[cpg] == 1:
            tissue = random.choice(tissues)
            tissue_dict[tissue].append(cpg)

    return tissue_dict


def create_reference_array(tissue_dict, cpg_dict, mean, sd):
    ref_array_lists = [[] for i in range(len(tissue_dict))]

    for cpg in range(len(cpg_dict)):
        print(cpg)
        if cpg_dict[cpg] == 1:
            ref_array_lists[cpg] = np.random.normal(mean, sd, len(cpg_dict))
        else:
            l

        print(ref_array_lists)

if __name__ == "__main__":

    # accepts user input
    # all choices must be greater than 1
    # default is 5 individuals, 3 tissues, and 3 CpGs
    parser = argparse.ArgumentParser()
    parser.add_argument("--individuals", help="number of individuals", type=int, choices=range(2,100), default=5)
    parser.add_argument("--tissues", help="number of tissues", type=int, choices=range(2,100), default=3)
    parser.add_argument("--cpgs", help="number of CpG sites", type=int, choices=range(2,100), default=3)

    args = parser.parse_args()

    individual = args.individuals
    cpg = args.cpgs
    tissue = args.tissues

    # hardcoded mean and sd from methylation dataset mentioned in refactor paper
    # TODO change this to something better
    mean = 7000
    sd = 435

    # TODO make this into function, not in main method
    cpg_dict = choose_dmr(cpg)
    tissue_dict = assign_cpg_to_tissue(cpg_dict, tissue)
    create_reference_array(tissue_dict, cpg_dict, mean, sd)