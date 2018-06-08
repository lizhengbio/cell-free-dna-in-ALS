import csv
import subprocess
import pysam
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def get_sequence(chrom, start, end, range):
    """
    using a compressed hg19 file, find the sequence for a region of interest
    :param chrom: chromosome of site of interest
    :param start: range start
    :param end: range end
    :return: sequence
    """
    # print(start)
    # print(end)
    # start = int(start) - range 
    # end = int(end) + range
    # print(start)
    # print(end)

    # command to get the region from the two bit file from fasta
    cmd = ["/ye/zaitlenlabstore/christacaggiano/twoBit/twoBitToFa", "/ye/zaitlenlabstore/christacaggiano/twoBit/hg38.2bit",
           "stdout", "-seq=" + chrom, "-start=" + str(start), "-end=" + str(end)]

    # call command and get output
    result = subprocess.check_output(cmd)

    return result.decode().upper()


def format_fasta_entry(entry, cpg_dict, tissue, entry_type):

    header = ""

    for i in range(len(entry)):
        if entry[i] is not "\n":
            header += entry[i]
        else:
            break

    sequence = entry[i:].replace("\n", "")
    sequence = mark_cpgs(sequence, cpg_dict)

    header = header + "|" + "total_cpgs=" + str(sequence.count("*")) + "|" + tissue + "|" + entry_type 


    return header, sequence

def mark_cpgs(sequence, cpg_dict):

    offset = -1 
    for i in cpg_dict: 
        location = i + offset
        if sequence[location] == "C":
            sequence = sequence[:location+1] + "*" + sequence[location+1:]
            offset += 1

    return sequence

def convert_to_methylated(sequence):

    sequence = sequence.replace("*", "")

    for i in range(len(sequence)-1):
        if sequence[i] == "C" and sequence[i+1] != "G": 
            sequence = sequence[:i] + "T" + sequence[i+1:]
    return sequence

def count_region_cpgs(bed_file):

    cpg_count_dict = {}


    # for each cpg in the bed file of interest
    with open(bed_file) as bed:

        # next(bed)  # skip header

        for line in bed:

            line = line.strip("\n")
            chrom, start, end = line.split("\t")[0], line.split("\t")[1], line.split("\t")[2]  # get pos information
            cpg_start = line.split("\t")[4]

            cpg_distance = int(cpg_start) - int(start)

            if (chrom, start, end) not in cpg_count_dict:
                cpg_count_dict[(chrom, start, end)] = [cpg_distance]
            else:
                cpg_count_dict[(chrom, start, end)].append(cpg_distance)

    return cpg_count_dict

def open_sam(sam_file):
    
    return pysam.AlignmentFile(sam_file, "r")

def pairwise_alignment(reference, read): 
    alignments = pairwise2.align.localms(reference, read,  5, -1, -10, -10)
    return(format_alignment(*alignments[0]))

    
if __name__ == "__main__":

    # TODO: command line arguments

    bed_file = "brain_hg38_125_results.txt"  # file that primers are being gotten for
    tissue = "brain"
    output_file = bed_file.replace(".txt", "_probes.txt")  # output file name
    sam_file = "/ye/zaitlenlabstore/christacaggiano/cfDNA_project/primary_methylation_data_ALS/1_BC_S1_L001_R1_001/unsortedButMerged_ForBismark_file/sorted.bam.bam"
    sam_file = open_sam(sam_file)

    cpg_count_dict = count_region_cpgs(bed_file)

   
    with open(output_file, "w") as f:

        for cpg in cpg_count_dict:
            if len(cpg_count_dict[cpg]) > 3:

                chrom, start, end = cpg
                sequence = get_sequence(chrom, start, end, 125)
                header, sequence = format_fasta_entry(sequence, cpg_count_dict[cpg], tissue, "reference")
                
                print(header, file=f)
                print(sequence, file=f)

                print(header.replace("reference", "methylated"), file=f)
                print(convert_to_methylated(sequence), file=f)

                reads = sam_file.fetch(region=chrom + ":" + start + "-" + end)
                read_num = 1 
                for read in reads:
                    print(header.replace("reference", "cfDNA_read_number=" + str(read_num)), file=f)
                    print(pairwise_alignment(convert_to_methylated(sequence), read.seq), file=f)
                    read_num+=1
                    

                 













