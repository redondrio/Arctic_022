import sys

file = sys.argv[1]

with open(file,'r') as infile, open(file.split('_pre')[0]+'.txt','w') as outfile:
	for line in infile:
		fig = line.split('\t')[0]
		domains = line.split('\t')[1].split(';')
		for domain in domains:
			code = domain.split(' ')[0]
			outfile.write('{}\t{}\n'.format(fig,code))



# Esto es para write (overwrite)
with open('file','w') as outfile:
	outfile.write('{}string'.format(values))

# Esto es para append
with open('file','a') as outfile:
	outfile.write('{}string'.format(values))

# Esto es para leer y escribir a la vez
with open('file','r') as infile, open('out','w') as outfile:
	for line in infile
		outfile.write('{}'.format(line))