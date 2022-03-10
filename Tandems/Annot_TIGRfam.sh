#####################################
# HMMER CONTRA PFAM+TIGRFAM
# TIGR04056 for SusC
#####################################

# Concatenar Pfam-A y TIGRfam
# Coger solo el TIGR04056
# Importante que NO incluya el // al final
# pero si al principio para separar
tail -n 457642 tigrfam.hmm | head -2970 > tigr_susc_only.hmm
cd /media/disk2/aredondo/artico022_DNA/rds/fastas
cat /media/disk7/fer/SqueezeMeta/db/Pfam-A.hmm tigr_susc_only.hmm > Pfam-A_and_TIGRfam.hmm

# ESTO NO ESTÁ HECHO, SOLO LO PLANTEO
'''
# Para los TIGR04183 y TIGR04131 (PF18962)
# Estas son las señales de secreción del T9SS
# que quiero buscar en las CAZymes

tail -n 359785 tigrfam.hmm | head -243 > tigr_T9SS_A.hmm #TIGR04183 incluye //
tail -n 397408 tigrfam.hmm | head -282 > tigr_T9SS_B.hmm #TIGR04131 incluye //

cat tigr_T9SS_A.hmm tigr_T9SS_B.hmm > tigr_T9SS.hmm
cp Pfam-A_and_TIGRfam.hmm Pfam-A_and_TIGRfam2.hmm
cat /media/disk2/aredondo/db/Pfam-A_and_TIGRfam2.hmm tigr_T9SS.hmm > Pfam-A_and_TIGRfam.hmm
rm Pfam-A_and_TIGRfam2.hmm
rm tigr_T9SS.hmm
'''

# Copia de seguridad del conf.pl
cp SqueezeMeta_conf.pl original_SqueezeMeta_conf.pl
# Poner el path en SqueezeMeta_conf.pl
$pfam_db   = "/media/disk2/aredondo/artico022_DNA/rds/fastas/Pfam-A_and_TIGRfam.hmm"
# Activar nokegg y nocog para que el script 7 solo reanote PFAM
$nocog           = 1
$nokegg          = 1

# Renombrar los outputs
cd intermediate
mv 05.artico022_DNA.pfam.hmm original_05.artico022_DNA.pfam.hmm 
cd results
mv 07.artico022_DNA.fun3.pfam original_07.artico022_DNA.fun3.pfam

# Correr HMMer
# Hay que hacerlo desde el directorio superior al proyecto
cd ..
conda activate SQM
nohup 05.run_hmmer.pl artico022_DNA &   # run hmmer
nohup 07.fun3assign.pl artico022_DNA &  # parse hmmer output
# recuerda renombrar los 07. que no son pfam después
# de correr el script, porque si no SQMtools se lía
# hay que sustituir el fun3 

# These steps are needed for SQMtools
nohup 12.funcover.pl artico022_DNA &    # calculate coverages
nohup 13.mergeannot2.pl artico022_DNA & # merge for SQMtools

#####################################
# EXTRAER ANOTACIONES DE TIGRFAM
# TIGR04131/83 for T9SS signal
#####################################

# Get the orfs that have TIGR04131/83 annotated
grep 'TIGR041[31|83]' 07.artico022_DNA.fun3.pfam > t9ss.tigr

# Field separator \t
# While in the first file (FNR==NR),
# load the line ($0) in a dictionary with the orf ID ($1) as key
# Once in the second file,
# If the orf was not in the first file (a[$1]=="")),
# keep the line in the old file ($0), else
# print the line from the new file (a[$1])

awk -F'\t' 'FNR==NR{a[$1]=$0;next;} {if (a[$1]=="") print $0; else print a[$1]}' <(cat 07.artico022_DNA.fun3_plus.tigr) <(cat 07.artico022_DNA.fun3.pfam) > 07.artico022_DNA.fun3.pfamnew

# These steps are needed for SQMtools
nohup 12.funcover.pl artico022_DNA &    # calculate coverages
nohup 13.mergeannot2.pl artico022_DNA & # merge for SQMtools

















