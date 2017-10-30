# 29 Sept
# python script to call WGBS processes


import argparse
import subprocess
import os


if __name__=="__main__":

    # takes in command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("run", help="path containing fastq files")
    parser.add_argument("output", help="path where outputdir should be created")
    parser.add_argument("file", help="fastq file for processing")
    args = parser.parse_args()
    output = args.output
    run = args.run
    file = args.file

    if not(os.path.exists(output)):
        os.mkdir(output)

    # split fastq files into groups of 18 million reads
    subprocess.call(["fastq_split.sh", file, output], shell=True)





