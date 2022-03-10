# Mincount for CAZy annots
# This receives a file created with the CAZyDB_ref_annots.sh script
# and outputs a file with the annotations that have at least 'threshold' counts

import sys

file = sys.argv[1]
threshold = int(sys.argv[2])

with open(file,'r') as infile, open(file.split('.uniq')[0]+'.min'+str(threshold)+file.split('.uniq')[1],'w') as outfile:
	for line in infile:
		count = int(line.split(' ')[1])
		annot = line.split(' ')[2].strip('\n')
		if count >= threshold: outfile.write(annot+'\n')