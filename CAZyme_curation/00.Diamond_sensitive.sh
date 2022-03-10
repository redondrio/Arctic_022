# Hmmer_dbCAN
wget http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V9.txt
hmmpress dbCAN-HMMdb-V9.txt

mv artico022_DNA/results/05.artico022_DNA.pfam.hmm artico022_DNA/results/05.artico022_DNA.pfam_susc.hmm
less artico022_DNA/SqueezeMeta_conf.pl #changed the Pfam link to dbCAN-HMMdb-V9.txt
conda activate SQM1_4
nohup 05.run_hmmmer.pl artico022_DNA &

# La idea después es coger las secuencias que dbCAN anote y hacerle el diamond sensitive
# Habría que extraerlas del 03.artico022_DNA.faa

# Diamond_sensitive
# Se usan los parámetros del grupo de Bremen
nohup /home/aredondo/miniconda3/envs/SQM1_4/SqueezeMeta/bin/diamond blastp \
--sensitive -q results/03.artico022_DNA.faa -p 22 -d /media/disk2/aredondo/db/CAZy.dmnd \
-e 0.00000000000000000001 --id 30 --query-cover 40 --quiet -b NF \
-f 6 qseqid qlen sseqid slen pident length qcovhsp evalue bitscore qstart qend sstart send \
-o intermediate/04.artico022_DNA.CAZy_filtered.diamond &