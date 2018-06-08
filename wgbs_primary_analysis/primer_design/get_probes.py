import csv
import subprocess


def get_sequence(chrom, start, end):
    """
    using a compressed hg19 file, find the sequence for a region of interest
    :param chrom: chromosome of site of interest
    :param start: range start
    :param end: range end
    :return: sequence
    """

    start = int(start) - 250  # hardcoded range around the cpg start and end :(
    end = int(end) + 250

    # command to get the region from the two bit file from fasta
    cmd = ["/Users/Christa.Caggiano/Documents/software/twoBitToFa", "/Users/Christa.Caggiano/Documents/software/hg19.2bit",
           "stdout", "-seq=" + chrom, "-start=" + str(start), "-end=" + str(end)]

    # call command and get output
    result = subprocess.check_output(cmd)

    return result.decode().upper()


def format_fasta_entry(entry, cpg_number, tissue):

    header = ""

    for i in range(len(entry)):
        if entry[i] is not "\n":
            header += entry[i]
        else:
            break

    header = header + "|" + "total_cpgs=" + str(cpg_number) + "|" + tissue

    sequence = entry[i:].replace("\n", "")

    return header, sequence




if __name__ == "__main__":

    # TODO: command line arguments

    bed_file = "brain/brain_cpgs_50_results.txt"  # file that primers are being gotten for
    tissue = "brain"
    output_file = bed_file.replace(".txt", "_probes.txt")  # output file name
    cpg_count_dict = {}

    # for each cpg in the bed file of interest
    with open(bed_file) as bed:

        next(bed)  # skip header

        for line in bed:

            line = line.strip("\n")
            chrom, start, end = line.split("\t")[0], line.split("\t")[1], line.split("\t")[2]  # get pos information

            previous_chrom = chrom
            previous_start = start

            if (chrom, start, end) not in cpg_count_dict:
                cpg_count_dict[(chrom, start, end)] = 1
            else:
                cpg_count_dict[(chrom, start, end)] += 1

    with open(output_file, "w") as f:

        for cpg in cpg_count_dict:
            if cpg_count_dict[cpg] > 2:

                chrom, start, end = cpg
                sequence = get_sequence(chrom, start, end)
                header, sequence = format_fasta_entry(sequence, cpg_count_dict[cpg], tissue)

                print(header, file=f)
                print(sequence, file=f)
