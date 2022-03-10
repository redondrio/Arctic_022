# Curation of IDs
# This receives a file created with the CAZyDB_ref_annots.sh script
# and outputs the selected annotations

import sys

file = sys.argv[1]

with open(file,'r') as infile, open(file.split('.uniq')[0]+'.curated'+file.split('.uniq')[1],'w') as outfile:
	for line in infile:
		print(line)
		accept = input("Accept? (Y/END/UNDO) ")
		if accept == "UNDO" : 
			outfile.write(memory+'\n')
			accept = input("Accept? (Y/END) ")
		if accept == "END" : break
		if accept == "Y" : outfile.write(line+'\n')
		memory = line