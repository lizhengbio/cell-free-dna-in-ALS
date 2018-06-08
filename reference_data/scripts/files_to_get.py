toget = open("files (1).txt", "r")
current = open("currentlyhave.txt", "r")

# get_list = []
# for i in current: 
# 	for j in current: 
# 		if i not in j: 
# 			get_list.append()

toget_list = []
for line in toget:
	if "fastq.gz" in line:
		toget_list.append(line[59:])

print(toget_list)
print(have_list)

s1 = set(toget_list)
s2 = set(have_list)

diff = s1.difference(s2)

file_names = [] 
for item in diff: 
	string_to_write = "https://www.encodeproject.org/files/" + item[:11] + "/" +"@@download/" + item 
	file_names.append(string_to_write)


with open("toget_list.txt", 'w') as f:
    for item in file_names:
    	f.write(item)
