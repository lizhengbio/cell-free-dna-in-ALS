# 22 Sept 17 
# author <christa.caggiano@ucsf.edu> 
# generates file list for wget -i mass data download 

files = []
url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP051/SRP051242/SRR"

def generate_urls(files, url, start, stop):
	for i in range(start, stop):
		files.append(url + str(i) + "/SRR" + str(i) + ".sra") 

generate_urls(files, url,  1720589, 1720752)
generate_urls(files, url,  2533657, 2533675 )

print(len(files))

with open("files.txt", 'w') as f:
	for url in files: 
		f.writelines(url)
		f.writelines('\n')
f.close()