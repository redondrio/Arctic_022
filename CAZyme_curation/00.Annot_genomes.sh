#####################################
# ANOTACIÓN DE GENOMAS INDIVIDUALES
#####################################

# Esto está metido en un script que se llama
# annotate_genome.sh
# Se ejecuta directamente como 
# nohup bash annotate_genome.sh &

for genomecode in $(cut -d ' ' -f 1 genome_list.txt)
do
	# ir al padre del directorio con los archivos
	cd /media/disk2/aredondo/artico022_DNA/genomes
	# descargar .faa
	wget ftp://ftp.patricbrc.org/genomes/$genomecode/$genomecode.PATRIC.faa
	# crear file.samples cuyo nombre es el del proyecto
	echo $genomecode$'\t'$genomecode.PATRIC.faa$'\tpair1' > $genomecode.samples
	# crear un proyecto vacío
	# para el proyecto vacío vale un archivo existente cualquiera
	SqueezeMeta.pl -m sequential -s $genomecode.samples -f . -empty
	# colocar el .faa en resultados para las anotaciones
	# cambiar el nombre para que SQM lo reconozca
	mv $genomecode.PATRIC.faa $genomecode/results/03.$genomecode.faa
	# el fichero SqueezeMeta_conf.pl define las variables del proyecto
	# al final, en la sección de opciones, define la variable $opt_db
	# esa variable es el directorio al mydblist.txt con los directorios de las DBs
	echo '$opt_db             = "/media/disk2/aredondo/artico022_DNA/genomes/mydblist.txt";' >> $genomecode/SqueezeMeta_conf.pl
	# ahora se ejecutan las anotaciones funcionales
	# como las bases de datos están actualizadas, hay
	# que usar el diamond de la carpeta /opt
	04.rundiamond.pl $genomecode
	05.run_hmmer.pl $genomecode
	07.fun3assign.pl $genomecode
done

cat genome_list.txt >> .done_genomes.txt
echo "" > genome_list.txt

# PARA AÑADIR MÁS BASES DE DATOS PERSONALIZADAS
# El directorio de la base de datos hay que meterlo en el fichero mydblist.txt
# La base de datos tiene que estar en formato .dmnd
## echo CAZyDB.dmnd$'\t'/media/disk2/aredondo/artico022_DNA/rds/CAZyDB.dmnd > mydblist.txt

# Dejar algo corriendo y salir
# Ten en cuenta que en el archivo genome_list tienen que estar los códigos de genomas NO analizados
#nohup bash annotate_genome.sh &
