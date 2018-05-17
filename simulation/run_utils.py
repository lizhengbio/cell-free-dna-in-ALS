import pandas as pd


def generate_matrices(intersect_file):
    data = pd.read_table(intersect_file, low_memory=False)

    meth_array = data["meth"].values
    unmeth_array = data["unmeth"].values
    reference = data.drop(["chr", "start", "end", "cpg", "chr.1", "start.1", "end.1", "meth", "unmeth", "total"], axis=1)

    return reference


def read_tissues(tissue_file):
    tissues = {}
    with open(tissue_file) as f:
        for line in f:
            line = line.strip("\n").split("\t")
            tissues[line[0]] = [x for x in line[1:] if not x == ""]
    return tissues


def calculate_tissue_means(reference,  tissue_dict):
    for tissue in tissue_dict:
        df = reference[tissue_dict[tissue]]
        print(type(df))
        print(df.mean(axis=0, numeric_only=True))


reference = generate_matrices("individual_1/intersect.txt")
tissues = read_tissues("tissues.txt")
calculate_tissue_means(reference, tissues)