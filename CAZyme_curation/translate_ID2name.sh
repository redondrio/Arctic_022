# Translate COG
# Just set file to the file to translate
# Make sure that the correct field is extracted

old_IFS=$IFS
IFS=$'\n'
for line in $(<$file)
do 
	id=$(echo $line | cut -d ' ' -f 3) #check this is correct for the file
	trans=$(grep $id /home/aredondo/Arctic/genomes/coglist.txt)
	echo $line $trans " " >> "trans_$file"
done
IFS=$old_IFS

# Translate KEGG
# Just set file to the file to translate
# Make sure that the correct field is extracted

old_IFS=$IFS
IFS=$'\n'
for line in $(<$file)
do 
	id=$(echo $line | cut -d ' ' -f 3) #check this is correct for the file
	trans=$(grep $id /home/aredondo/Arctic/genomes/keggfun2.txt)
	echo $line $trans " " >> "trans_$file"
done
IFS=$old_IFS

