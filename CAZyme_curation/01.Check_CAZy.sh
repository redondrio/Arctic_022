######################################
# CHECK CAZY ANNOTATIONS IN GENOMES
# AND EXTRACT TP, FP AND FN
######################################

#######################
# FIRST PART
#######################
# This script will extract CAZy annotations from the reference genomes
# and extract the CAZy annotations that SqueezeMeta did
# with the IDs that the PATRIC genome has for proteins
# Then, it will compare both lists of IDs and extract true positives,
# false positives and false negatives, taking the reference annotations
# as the real positivies

# How to run
# for id in $(cat .done_genomes.txt | cut -f 1); do bash CAZyme_annot_check.sh "GT" $id; done
# This script will focus on GTs
cazyfam=$1

# This script is tought to be run in the Arctic/genomes directory
# iterating over the different genome directories
# For now it is only done with one genome
genomeid=$2

echo "Begin" $genomeid "genome"
cd $genomeid/results
mkdir -p "$cazyfam"_quality

# File cazome.txt contains the annotated seqs in CAZy
# This will create file ref_cazyfam_ids.txt with the unique IDs of the reference CAZymes
# and remove it in case is was previously created

rm "$cazyfam"_quality/ref_"$cazyfam"_ids.txt

# Get the genes annotated as the family in CAZy 
grep "$cazyfam" cazome.txt | cut -f 1 > temp_ids.txt

old_IFS=$IFS
IFS=$'\n'

for id in $(<temp_ids.txt)
do
	echo $id
	grep -w "$id" 03.$genomeid.faa | cut -d ' ' -f 1 | tr -d ">" >> temp_fig.txt
done

IFS=$old_IFS

cat temp_fig.txt | sort > "$cazyfam"_quality/ref_"$cazyfam"_ids.txt

rm temp_*.txt
cd "$cazyfam"_quality

# File 07.$genomeid.fun3.CAZyDB.dmnd contains the SQM annotations
# This will create file sqm_"$cazyfam"_ids.txt with the unique IDs for the annotated CAZymes
grep "$cazyfam" ../07.$genomeid.fun3.CAZy > sqm_"$cazyfam".txt
cut -f 1 sqm_"$cazyfam".txt | sort | uniq > sqm_"$cazyfam"_ids.txt

# File 04.313590.8.CAZyDB.dmnd.diamond contains the diamond output
# This will create file dmnd_"$cazyfam".tsv with the diamond-BlastX results with all the quality scores
# This file will be loaded in R to see how different thresholds change the accuracy
grep "$cazyfam" ../../intermediate/04.$genomeid.CAZy.diamond > dmnd_"$cazyfam".tsv

# Check TP, FP and FN comparing the reference IDs and the annotated IDs
# Extract IDs common to both, considered true positives
comm -12 ref_"$cazyfam"_ids.txt sqm_"$cazyfam"_ids.txt > TP_"$cazyfam".txt
# Extract IDs only annotated by SQM, considered false positives
comm -13 ref_"$cazyfam"_ids.txt sqm_"$cazyfam"_ids.txt > FP_"$cazyfam".txt
# Extract IDs only annotated in the reference, considered false negatives
comm -23 ref_"$cazyfam"_ids.txt sqm_"$cazyfam"_ids.txt > FN_"$cazyfam".txt

#######################
# SECOND PART
#######################
# This part will check the complete annotations of FP and FN proteins,
# to see if the unexpected result is due to a missannotation of the reference
# or an error of SqueezeMeta.
# For that, KEGG, COG and PFAM annotations will be extracted for all FP and FN

# Extract the rest of annotations by other databases
# the -w option is to avoid fig.3 to match fig.3 and fig.30

# Because the loop needs to append, files must be deleted first
# if the script is to be run twice
rm ??_"$cazyfam"_*.txt
echo "Extracting associated annotations"
for fig in $(cat FP_"$cazyfam".txt)
do
	grep -w $fig ../07.$genomeid.fun3.cog >> FP_"$cazyfam"_cog.txt
	grep -w $fig ../07.$genomeid.fun3.pfam >> FP_"$cazyfam"_pfam_pre.txt
	grep -w $fig ../07.$genomeid.fun3.kegg >> FP_"$cazyfam"_kegg.txt
done
for fig in $(cat FN_"$cazyfam".txt)
do
	grep -w $fig ../07.$genomeid.fun3.cog >> FN_"$cazyfam"_cog.txt
	grep -w $fig ../07.$genomeid.fun3.pfam >> FN_"$cazyfam"_pfam_pre.txt
	grep -w $fig ../07.$genomeid.fun3.kegg >> FN_"$cazyfam"_kegg.txt
done
for fig in $(cat TP_"$cazyfam".txt)
do
	grep -w $fig ../07.$genomeid.fun3.cog >> TP_"$cazyfam"_cog.txt
	grep -w $fig ../07.$genomeid.fun3.pfam >> TP_"$cazyfam"_pfam_pre.txt
	grep -w $fig ../07.$genomeid.fun3.kegg >> TP_"$cazyfam"_kegg.txt
done

# Translate KEGG, COG and PFAM IDs to names
# The translation files are already in the Arctic directory
# and they are keggfun2.txt, coglist.txt y pfam.dat
# The idea is, for every line, add the ID and the translation

# This change is for the files to be read by line
old_IFS=$IFS
IFS=$'\n'

# KEGGs
echo "Translating"
for file in $(ls ??_"$cazyfam"_kegg.txt)
do
for line in $(<$file)
do 
	id=$(echo $line | cut -f 1)
	trans=$(grep $(echo $line | cut -f 2) ../../../keggfun2.txt)
	echo $id $trans >> $(basename -s .txt $file)_names.txt
done
done

# COGs
for file in $(ls ??_"$cazyfam"_cog.txt)
do
for line in $(<$file)
do 
	id=$(echo $line | cut -f 1)
	trans=$(grep $(echo $line | cut -f 2) ../../../coglist.txt)
	echo $id $trans >> $(basename -s .txt $file)_names.txt
done
done

IFS=$old_IFS

# Check CAZY annotation of FP
# These have been annotated as CAZymes, but should have not
# Let's check how their annotation was, maybe it was really
# really good, and then the reference is the one missing

for fig in $(cat FP_"$cazyfam".txt)
do
	grep -w $fig ../03.$genomeid.faa >> FP_"$cazyfam"_reference.txt
done

# Check CAZY annotation of FN
# These have not been annotated as CAZymes, but should have
# Let's check how their annotation was, maybe it was bad and
# then the reference is the one missing

for fig in $(cat FN_"$cazyfam".txt)
do
	grep -w $fig ../03.$genomeid.faa >> FN_"$cazyfam"_reference.txt
done

# Format PFAM annotations
for file in $(ls ??_"$cazyfam"_pfam_pre.txt)
do
	python3 ../../../pfam_parsing.py $file
done

echo "Genome" $genomeid "done"
cd ../../..

















