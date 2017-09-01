# ramachandran Plots summer boot camp project
# 29 August 2017
# author <christa.caggiano@ucsf.edu>

import argparse
import math
import matplotlib.pyplot as plt


######## vector operations #############

def add_vector(vector1, vector2):
    """
    add two vectors
    :param vector1:
    :param vector2:
    :return: a vector
    """
    return vector1[0] + vector2[0], vector1[1] + vector2[1], vector1[2] + vector2[2]


def subtract_vector(vector1, vector2):
    """
    subtract two vectors
    :param vector1:
    :param vector2:
    :return: a vector
    """
    return (vector1[0] - vector2[0]), (vector1[1] - vector2[1]), (vector1[2] - vector2[2])


def normalize(vector):
    """
    normalize a vector in the direction of the unit vector
    :param vector:
    :return: a normalize vector
    """
    sum = 0
    for value in vector:
        sum += value**2
    return [v/math.sqrt(sum) for v in vector]


def dot_product(vector1, vector2):
    """
    dot product of two vectors
    :param vector1:
    :param vector2:
    :return: an int
    """
    return vector1[0] *vector2[0] + vector1[1]*vector2[1] + vector1[2] * vector2[2]


def cross_product(vector1, vector2):
    """
    take the cross product of a 3-d vector
    :param vector1:
    :param vector2:
    :return: a vector
    """
    return vector1[1]*vector2[2] - vector1[2]*vector2[1], vector1[2]*vector2[0] - vector1[0]*vector2[2], \
           vector1[0]*vector2[1] - vector1[1]*vector2[0]


######### coordinate extraction ############

def get_coordinates(pdb_file, chain):
    """
    extracts the cartesian coordinates from a pbd file
    :param pdb_file: read in file containing pdb information
    :param chain: which chain of the protein the user wants plotted
    :return: dictionary of all coordinates for a particular atom name
    """
    keys = ["C", "N", "CA"]  # atom names of interest
    coordinate_dict = {key:[] for key in keys}
    read_file = pdb_file.read().splitlines()  # split the line in the pdb file to remove all whitespace

    for i, line in enumerate(read_file):
        entry = line.split()

        # if the entry is an ATOM entry for the appropriate chain, and the molecule is of interest add to dict
        if entry[0] == "ATOM" and entry[4] == chain and entry[2] in coordinate_dict:
            coordinate_dict[entry[2]].append((float(entry[6]), float(entry[7]), float(entry[8])))

    return coordinate_dict


def make_vectors(key1, key2, coord_dict, i1, i2):
    """
    makes vector from cartesian coordinates through vector subtraction)
    :param key1: which molecule should be chosen
    :param key2:
    :param coord_dict: dict holding all the cartesian coordinates
    :param i1: fetches the appropriate entry from the dict
    :param i2:
    :return: a vector
    """
    return subtract_vector(coord_dict[key1][i1], coord_dict[key2][i2])

######### calculate torsion angles ############

def phi(coord_dict, i):
    """
    extracts the coordinates for a phi calculation
    :param coord_dict:
    :param i:
    :return: a tuple of three vectors
    """
    return make_vectors("C", "N", coord_dict, i-1, i), make_vectors("N", "CA", coord_dict, i, i), make_vectors("CA", "C", coord_dict, i, i)


def psi(coord_dict, i):
    """
    extracts appropriate coordinates for a psi calculation
    :param coord_dict:
    :param i:
    :return: a tuple of three vectors
    """
    return make_vectors("N", "CA", coord_dict, i, i), make_vectors("CA", "C", coord_dict, i, i), make_vectors("C", "N", coord_dict, i, i+1)


def return_angle(vector_group):
    """
    calculates the cross product between vectors of interest, takes the dot product, and uses the ratio between those
    dot products to calculate the torsion angle
    from rama algorithm included @
    :param vector_group:
    :return: angle in degrees
    """
    n1, n2 = normalize(cross_product(vector_group[0], vector_group[1])), normalize(cross_product(normalize(vector_group[1]), vector_group[2]))
    x = dot_product(n1, n2)
    y = dot_product(cross_product(n1, normalize(vector_group[1])), n2)
    return math.degrees(math.atan2(y, x))


# TODO make this handle exceptions
def calculate_phi_psi(coord_dict, angle_type):
    """
    controls calculations for either phi or psi given input to the fuction call
    :param coord_dict: dict of cartesian coordinates
    :param angle_type: either phi or psi
    :return: a list of torsion angles
    """
    torsion_list = []
    for i in range(1,len(coord_dict["C"])-1):
        if angle_type == "phi":
            torsion_list.append(return_angle(phi(coord_dict, i)))
        else:
            torsion_list.append(return_angle(psi(coord_dict, i)))
    return torsion_list


if __name__ == "__main__":

    # takes user input to the file
    parser = argparse.ArgumentParser()
    parser.add_argument("file_name", type=str, default="4jzy.pdb")  # file name must be included
    parser.add_argument("chain", type=str, default="A")  # chain information must be included

    args = parser.parse_args()

    # opens file without storing in memory
    with open(args.file_name, "r") as pdb_file:
        coordinate_dict = get_coordinates(pdb_file, args.chain)
        phi_list = calculate_phi_psi(coordinate_dict, "phi")
        psi_list = calculate_phi_psi(coordinate_dict, "psi")

        # makes rama plot
        plt.ylim(-180, 180)
        plt.xlim(-180, 180)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.grid()
        plt.plot(phi_list, psi_list, 'go')
        plt.show()
