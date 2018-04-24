# gets methylation specific primers meeting input criteria
# primer analysis from MethPrimer2 http://www.urogene.org/cgi-bin/methprimer2/MethPrimer.cgi
# March 2018
# author <christa.caggiano@ucsf.edu>

# imports - requires several programs for interacting with webpages
import requests
import subprocess
import csv
from bs4 import BeautifulSoup


def has_primer(sequence):
    """
    checks whether the sequence has a primer using MethPrimer2
    :param sequence: string of nucleotides
    :return: primer status, output information
    """

    # MethPrimer2 parameters
    payload = {
        "SEQUENCE": sequence,
        "PRIMER_TASK": "PICK_MSP_PRIMERS",
        "CPG_ISLAND_SIZE": 100,
        "CPG_SHIFT": 1,
        "CPG_OE": 0.6,
        "CPG_ISLAND_PERCENT": 50.0,
        "BLAT_GENOME": "",
        "PRIMER_SEQUENCE_ID":  "",
        "PRIMER_NUM_RETURN": 5,
        "TARGET_REGION": "",
        "EXCLUDED_REGION": "",
        "PRIMER_MIN_TM": 59,
        "PRIMER_OPT_TM": 60,
        "PRIMER_MAX_TM": 61,
        "PRIMER_MIN_SIZE": 20,
        "PRIMER_OPT_SIZE": 25,
        "PRIMER_MAX_SIZE": 30,
        "PRIMER_MIN_GC": 10,
        "PRIMER_OPT_GC": 40,
        "PRIMER_MAX_GC": 80,
        "PRIMER_GC_CLAMP": 0,
        "PRIMER_CS": 4,
        "PRIMER_MAX_POLY_X": 5,
        "PRIMER_MAX_POLY_T": 8,
        "PRIMER_MAX_TM_DIFF": 5,
        "PRODUCT_MIN_SIZE": 50,
        "PRODUCT_OPT_SIZE": 150,
        "PRODUCT_MAX_SIZE": 180,
        "PRODUCT_CG_MIN": 4,
        "PRIMER_CG3_POSITION": 3,
        "PRIMER_CGS": 1,
        "SET_TM_DIFF": 2,
        "ENZYME_TYPE": "XhoI:CTCGAG",
        "ENZYME_MIN_CNT": 1,
        "ENZYME_MIN_TERINUS_DIST": 30,
        "PROBE_MIN_SIZE": 15,
        "PROBE_OPT_SIZE": 20,
        "PROBE_MAX_SIZE": 30,
        "PROBE_MIN_GC_PERCENT": 30,
        "PROBE_OPT_GC_PERCENT": 50,
        "PROBE_MAX_GC_PERCENT": 80,
        "PROBE_MAX_GC": 2,
        "PROBE_POLY_X": 5,
        "PROBE_HIGER_TM": 10,
        "PROBE_MIN_NON_CPG": 1,
        "PROBE_OPT_NON_CPG": 3,
        "PROBE_MAX_NON_CPG": 5,
        "PROBE_MIN_CPG": 1,
        "PROBE_OPT_CPG": 3,
        "PROBE_MAX_CPG": 5,
        "PRODUCT_NESTED_MIN_SIZE": 100,
        "PRODUCT_NESTED_OPT_SIZE": 300,
        "PRODUCT_NESTED_MAX_SIZE": 500,
    }

    # post request to MethPrimer2
    r = requests.post("http://www.urogene.org/cgi-bin/methprimer2/MethPrimer_result.cgi", data=payload)
    text = r.text

    # return empty set if there is no primer found
    if "No Acceptable Primers Were Found!" in text:
        return []

    else:

        soup = BeautifulSoup(r.content, "lxml")
        table = soup.find_all("table", class_="product")[0]  # get primer output

        important_rows = [1, 2, 6, 7]  # rows containing information about methylation specific primers
        important_items = [0, 1, 3, 4, 10]  # sequence information, tm information

        rows = list(table.find_all("tr"))  # get the primer info from website in an easier to use format
        results = []

        # for each row and column in the table of results, get important items
        for row_idx in important_rows:
            row = rows[row_idx]
            items = list(row.find_all("td"))

            for idx in important_items:
                results.append(items[idx].text)

        return results  # return the desired output to write out as a file


def get_sequence(chrom, start, end):
    """
    using a compressed hg19 file, find the sequence for a region of interest
    :param chrom: chromosome of site of interest
    :param start: range start
    :param end: range end
    :return: sequence
    """

    start = int(start) - 100  # hardcoded range around the cpg start and end :(
    end = int(end) + 100

    # command to get the region from the two bit file from fasta
    cmd = ["/Users/Christa.Caggiano/Documents/software/twoBitToFa", "/Users/Christa.Caggiano/Documents/software/hg19.2bit",
           "stdout", "-seq=" + chrom, "-start=" + str(start), "-end=" + str(end)]

    # call command and get output
    result = subprocess.check_output(cmd)

    return result


def format_result(sequence):
    """
    get sequence from two bit to fasta command which is a little wonky and not a nice string
    :param sequence: unformmated sequence
    :return: a sequence in a nice string with consistent casing
    """

    return sequence.decode().upper()


if __name__ == "__main__":

    # TODO: command line arguments

    bed_file = "pc-site.txt"  # file that primers are being gotten for

    output_file = bed_file.replace(".txt", "_primers.txt")  # output file name

    # for each cpg in the bed file of interest
    with open(bed_file) as bed, open(output_file, "w") as out:

        next(bed)  # skip header

        out_writer = csv.writer(out, delimiter="\t")

        for line in bed:

            line = line.strip("\n")
            chrom, start, end = line.split("\t")[0], line.split("\t")[1], line.split("\t")[2]  # get pos information
            sequence = get_sequence(chrom, start, end)  # get sequence
            results = has_primer(format_result(sequence))  # call MethPrimer2 and get primer

            if results:  # if there is a result, output it 
                out_writer.writerow([chrom, start, end] + line.split("\t")[3:] + results)
