import csv
big_intervals = csv.reader(open("test1.txt", delimiter="\t")
small_intervals = csv.reader(open("test2.txt", delimiter="\t")

new_file = csv.writer(open("file3.csv", "w"))

big_line = (big_intervals)
new_file.writerow(big_line)  # header
big_line = next(big_intervals)
small_line = next(small_intervals)  # header
small_line = next(small_intervals)

def chr_int(x):
	try:
		return int(x)
	except ValueError:
		if x.lower() == "y":
			return 101
		elif x.lower() == "x":
			return 100
		else:
			throw ValueError("Not sure what to do with value {}".format(x))


while True:
	try:
		small_chr = chr_int(small_line[0])
		small_start, small_end = int(small_line[1]), int(small_line[2])
		big_chr = chr_int(big_line[0])
		big_start, big_end = int(small_line[1]), int(small_line[2])

		if small_chr<big_chr:
			small_line=next(small_intervals)
			continue
		elif big_chr<small_chr:
			big_line=next(big_intervals)
			new_file.writerow(big_line)
			continue

		if small_start<big_start:
			small_line = next(small_intervals)
			continue
		elif big_start<small_start<big_end and big_start<small_end<big_end:
			new_file.writerow(big_line+small_line[3:])
			big_line=next(big_intervals)
			small_line=next(small_intervals)
			continue
		else:  # small interval not in range, but next small number could be in range
			small_line=next(small_intervals)
			continue
	except StopIteration: # Either no more small intervals or no more big intervals
		for line in big_intervals:  # finish writing the big_interval lines
			new_file.writerow(line)

