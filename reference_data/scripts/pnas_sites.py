# script that takes large methylatino data 
# and extracts sites from PNAS paper Lehman-Werner et al (2015)
# dec 2017 
# author <christa.caggiano@ucsf.edu> 
#@TODO make more robust to handle many file types

# imports 

# import argparse # not used in this version 


def check_sites(chr, start, pnas_sites):
    """ check if methylation file contains desired pnas site 
    @params: 
    chr: string of chr 
    start: int of site start 
    pnas_sites: dictionary of pnas_sites to be investigated 

    """

    # for each site in the pnas site dictionary check if our site 
    # is within the bounds. If so return true 
    for site in pnas_sites:
        if site[1] <= start <= site[2] and chr == site[0]:
            return True


def append_to_file(percent_meth):
    """ if a CpG site matches a PNAS site, append it to a file 
    @params: 
    percent_meth: CpG sites matching the PNAS sites  
    """
    f = open("pnas_sites.txt", "w") # write to file called pnas_sites

    for item in percent_meth:
        f.write("%s\n" % item) # write the percent for each item in our list 

if __name__ == "__main__":

    # # takes in command line arguments for input, output, and which file to use
    # # @TODO make run/output optional
    # parser = argparse.ArgumentParser()
    # parser.add_argument("input", help="path containing fastq files")
    # parser.add_argument("output", help="path where outputdir should be created")
    # parser.add_argument("file_name", type=int, help="fastq file for processing")
    # args = parser.parse_args()
    #

    # @TODO incorporate this information 
    # what tissue type is this 
    # sample_type = "cfDNA"
    # surrounding_sites = 6

    # files to check PNAS sites im 
    sample_list = ['results_2_trunc.txt']

    # pnas_sites to be checked * important hg18 sites 
    # @TODO take in a file so sites to b e checked can be flexible 
    pnas_sites = [("Chr18", 74692045, 7469221), ("Chr10", 79936919, 79937042), ("Chr10", 3283832, 3283996), ("Chr2", 79347448, 79347588)]

    # for each file to check, read in one line at a time and check whether it is a PNAS site  
    percent_meth = []
    for file in sample_list:
        with open(file) as f:
            
            for line in f:
                chr, start, meth, total = line.split()[0], int(line.split()[1]), int(line.split()[2]), int(line.split()[4]) # split line based on how I expect the file type to be formatted 
               
                # if the site exist, append it to our dictionary 
                if check_sites(chr, start, pnas_sites):
                    percent_meth.append((file, chr, start, meth/total))

    # append all sites to a file 
    append_to_file(percent_meth)





