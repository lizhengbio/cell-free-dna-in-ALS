from __future__ import print_function
import os

file_list = [line.strip() for line in open("file_list.txt", 'r')]
for file in file_list:
    with open(file) as f:
        os.system("rm -f " + file + "_cleaned.txt")
        with open(file + "_cleaned.txt", "a") as out:
            start_printing = False
            for line in f:
                if start_printing:
                    print(line.rstrip(), file=out)
                if line.startswith("\"ID_REF\""):
                    start_printing = True
                    print(line.rstrip(), file=out)
    print(open(file + "_cleaned.txt",).read())
