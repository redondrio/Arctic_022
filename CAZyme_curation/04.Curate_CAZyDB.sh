###############################
# GET THE COG/KEGG/PFAM IDS
# FROM CAZYDB ANNOTATIONS THAT
# WILL BE USED AS REFRENCE
###############################
# This part will het the COG/KEGG/PFAM annotations that were made on the whole CAZyDB
# to be used as reference annotations to identify other CAZymes.
# These will extracted in CAZyDB.uniq.* files.
# There is a script to filter by abundance (mincount.py),
# but the chosen method was manual curation (curation.py)

# COGs
# Remove comment lines
# Get only the BESTHIT column (consensus column not considered)
# Split semicolons, sort and get unique with counts and sort, most abundant first (remove extra blanks)

grep -v \# 07.CAZyDB.fun3.cog | grep '|GH' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.GH.uniq.cog
grep -v \# 07.CAZyDB.fun3.cog | grep '|GT' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.GT.uniq.cog
grep -v \# 07.CAZyDB.fun3.cog | grep '|PL' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.PL.uniq.cog
grep -v \# 07.CAZyDB.fun3.cog | grep '|CE' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.CE.uniq.cog
grep -v \# 07.CAZyDB.fun3.cog | grep '|CBM' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.CBM.uniq.cog

#wc -l CAZyDB.uniq.cog
#cut -d ' ' -f 2 CAZyDB.uniq.cog | paste -sd+ - | bc
# Got 4720 unique COG/ENOG IDs
# With 1641269 total counts

# KEGGs
# Remove comment lines
# Get only the BESTHIT column (consensus column not considered)
# Split semicolons, sort and get unique with counts and sort, most abundant first (remove extra blanks)

grep -v \# 07.CAZyDB.fun3.kegg | grep '|GH' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.GH.uniq.kegg 
grep -v \# 07.CAZyDB.fun3.kegg | grep '|GT' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.GT.uniq.kegg 
grep -v \# 07.CAZyDB.fun3.kegg | grep '|PL' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.PL.uniq.kegg
grep -v \# 07.CAZyDB.fun3.kegg | grep '|CE' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.CE.uniq.kegg
grep -v \# 07.CAZyDB.fun3.kegg | grep '|CBM' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.CBM.uniq.kegg

#wc -l CAZyDB.uniq.kegg
#cut -d ' ' -f 2 CAZyDB.uniq.kegg | paste -sd+ - | bc
# Got 1724 unique KEGG IDs
# With 1418405 total counts

# PFAM
# Remove comment lines
# Split semicolons, sort and get unique with counts and sort, most abundant first (remove extra blanks)

grep -v \# 07.CAZyDB.fun3.pfam | grep '|GH' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.GH.uniq.pfam 
grep -v \# 07.CAZyDB.fun3.pfam | grep '|GT' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.GT.uniq.pfam 
grep -v \# 07.CAZyDB.fun3.pfam | grep '|PL' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.PL.uniq.pfam
grep -v \# 07.CAZyDB.fun3.pfam | grep '|CE' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.CE.uniq.pfam
grep -v \# 07.CAZyDB.fun3.pfam | grep '|CBM' | cut -f 2 | tr ';' '\n' | sort | uniq -c | sort -nr -k 1 | tr -s ' ' > uniques/CAZyDB.CBM.uniq.pfam

#wc -l CAZyDB.uniq.pfam
#cut -d ' ' -f 2 CAZyDB.uniq.pfam | paste -sd+ - | bc
# Got 2221 unique PFAM IDs
# With 2481687 total counts

# This python script will extract only the annotations with at least 'threshold' counts
# and change the 'uniq' in the name for 'min100'
for file in $(ls *uniq*)
do
python3 mincount.py $file 500
done

# This python script will translate KEGG and COG for manual cuartion
# Translate COG
# Just set file to the file to translate
# Make sure that the correct field is extracted

old_IFS=$IFS
IFS=$'\n'
for file in $(ls *CBM.uniq.cog)
do
echo $file
for line in $(<$file)
do 
	id=$(echo $line | cut -d ' ' -f 3) #check this is correct for the file
	trans=$(grep $id /home/aredondo/Arctic/genomes/coglist.txt)
	echo $line $trans " " >> "trans_$file"
done
done

# Translate KEGG
# Just set file to the file to translate
# Make sure that the correct field is extracted

for file in $(ls *CBM.uniq.kegg)
do 
echo $file
for line in $(<$file)
do
	id=$(echo $line | cut -d ' ' -f 3) #check this is correct for the file
	trans=$(grep $id /home/aredondo/Arctic/genomes/keggfun2.txt)
	echo $line $trans " " >> "trans_$file"
done
done
IFS=$old_IFS

# Now do the manual curation
# for translated KEGG/COG
for file in $(ls trans_*)
do
	echo $file
	python3 curation.py $file
done

# for PFAM
for file in $(ls *uniq.pfam)
do
	echo $file
	python3 curation.py $file
done

# Now extract just the id
# cut -d ' ' -f 3 trans_* > *curid*











