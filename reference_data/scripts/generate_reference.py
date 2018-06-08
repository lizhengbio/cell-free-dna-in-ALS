# generate merged reference matrix from soft_matrix files obtained from GEO
# containing average beta values
# Jan 2018
# Christa Caggiano

import pandas as pd
from functools import reduce
from __future__ import print_function
import os
import argparse


def make_clean_file(file):

    with open(file) as f:

        clean_file_name = f.replace(".txt", "_cleaned.txt")
        remove_metadata(clean_file_name)


def remove_metadata(file):

    os.system("rm -f " clean_file)  # if a cleaned file exists, delete it and start fresh

    try:
        # open file and remove the metadata. Print to a new file.
        with open(clean_file, "w") as out:
            start_printing = False
            for line in f:
                if line.startswith("!Sample_geo_accession"):
                    print(line.rstrip(), file=out)
                elif start_printing:
                    print(line.rstrip(), file=out)
                elif line.startswith("\"ID_REF\""):
                    start_printing = True

    os.system("gzip " + file ) # zip up file with metadata because they tend to be large

    except Exception:

        # if the file is not in the proper format, throw an error
        print("ERROR: Please input a standard soft_matrix txt/tsv file.")
        print("File must include Sample_geo_accession and ID_REF lines")

        raise


def intiate_merge(file_list):
    df_list = ["df_" + str(i) for i in range(len(file_list))]

    for df in df_list:
        df = load_df(df)

    merge_df(df_list)


def load_df(df):
    try:
        return round_df(pd.read_table(df))
    except Exception:
        print("ERROR: unable to load data frame. Please make sure you have sufficient memory")
        raise


def round_df(df, decimals):
    return df.round(decimals=4)


def merge_df(df_list)

    df_final = reduce(lambda left, right: pd.merge(left, right, how='outer', on='!Sample_geo_accession'), df_list)
    df_final.to_csv('reference_matrix.txt', sep="\t")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("soft_file", type=str, help="txt file containing a list of soft matrix files to be merged")
    parser.add_argument(--non_standard, type=str, help="list of files not in soft matrix format


    soft_file_list = [line.strip() for line in open("soft_matrix_file_list.txt", 'r')]

    for file in file_list: make_clean_file(file)

    clean_file_list = [file + "_cleaned.txt" for file in file_list]
