#####################################
# BLAST de CAZyDB contra KEGG/COG/PFAM
#####################################

conda activate SQM

# ir al padre del directorio con los archivos
cd /media/disk2/aredondo/artico022_DNA/genomes
# traer el .faa a la carpeta
cp ../rds/CAZyDB.07312020.fa ./CAZyDB.faa
# crear file.samples cuyo nombre es el del proyecto
echo CAZyDB$'\t'03.CAZyDB.faa$'\tpair1' > CAZyDB.samples
# crear un proyecto vacío
# para el proyecto vacío vale un archivo existente cualquiera
SqueezeMeta.pl -m sequential -s CAZyDB.samples -f . -empty
# colocar el .faa en resultados para las anotaciones
mv 03.CAZyDB.faa CAZyDB/results/03.CAZyDB.faa
# ahora se ejecutan las anotaciones funcionales
nohup 04.rundiamond.pl CAZyDB &
nohup 05.run_hmmer.pl CAZyDB &
nohup 07.fun3assign.pl CAZyDB &