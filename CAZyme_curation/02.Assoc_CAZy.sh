######################################
# EXTRACT FUNCTIONS RELATED WITH
# CAZYMES IN THE GENOMES
######################################

# This extracts and counts all "$fun"s associated with CAZymes
# in any way (reference or SQM) to check what "$fun"s are more
# abundant in CAZymes

# This script is meant to be run from the genomes directory
# and receives the CAZy family to be annotated
# How to run: bash extract_assoc_fun.sh 'GT'
cazyfam=$1

for fun in "cog" "kegg" "pfam"
do
# Create directory and move there
mkdir assoc_"$cazyfam"_"$fun"
cd assoc_"$cazyfam"_"$fun"

# Extract all associated COG/KEGG annotations in the genomes
for id in $(cut ../.done_genomes.txt -f 1)
do
	cut -f 2 ../$id/results/"$cazyfam"_quality/FN_"$cazyfam"_"$fun".txt | sort >> FN_"$cazyfam"_"$fun"s.txt
done
for id in $(cut ../.done_genomes.txt -f 1)
do
	cut -f 2 ../$id/results/"$cazyfam"_quality/FP_"$cazyfam"_"$fun".txt | sort >> FP_"$cazyfam"_"$fun"s.txt
done
for id in $(cut ../.done_genomes.txt -f 1)
do
	cut -f 2 ../$id/results/"$cazyfam"_quality/TP_"$cazyfam"_"$fun".txt | sort >> TP_"$cazyfam"_"$fun"s.txt
done

# Count how many times annotation occurs
mkdir counts
for file in $(ls -p | grep -v /)
do
	cat $file | tr ';' '\n' | sort | uniq -c > counts/"count_$file"
done

# Sort them by descending frequency
mkdir sorted
for file in $(ls ./counts)
do
	cat ./counts/$file | sort -nr -k 1  | tr -s ' ' > sorted/"sorted_$file"
done

# Translate them to get the whole names
mkdir translated
old_IFS=$IFS
IFS=$'\n'
for file in $(ls ./sorted)
do
	echo $file
if [ $fun = "cog" ]
then
for line in $(<./sorted/$file)
do 
	id=$(echo $line | cut -d ' ' -f 3)
	trans=$(grep $id ../coglist.txt)
	echo $line $trans " " >> translated/"trans_$file"
done
elif [ $fun = "kegg" ]
then
for line in $(<./sorted/$file)
do 
	id=$(echo $line | cut -d ' ' -f 3)
	trans=$(grep $id ../keggfun2.txt)
	echo $line $trans "\n" >> translated/"trans_$file"
done
elif [ $fun = "pfam" ]
then
for line in $(<./sorted/$file)
do 
	id=$(echo $line | cut -d ' ' -f 3)
	trans=$(grep $id ../pfam.dat)
	echo $line $trans "\n" >> translated/"trans_$file"
done
fi
done
IFS=$old_IFS
# Go back to the /genomes directory
cd ..
done
