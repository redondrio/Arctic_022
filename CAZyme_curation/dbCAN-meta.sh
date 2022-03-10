#####
# INSTALACIÓN DE DBCAN-META
#####

wget https://bcb.unl.edu/dbCAN2/download/Tools/run_dbcan_09_16_2020.tar.gz
tar -xvf run_dbcan_09_16_2020.tar.gz run_dbcan/

# Instalo las cosas en el entorno SQM1_4
conda activate SQM1_4
pip install natsort

# Hay que preparar las bases de datos que se dan antes
diamond makedb --in CAZyDB.07312020.fa -d CAZy
# He usado la versión V9 de dbCAN-HMM pero he puesto V8 porque el script no estaba actualizado
hmmpress dbCAN-HMMdb-V8.txt
hmmpress stp.hmm
hmmpress tf-1.hmm
hmmpress tf-2.hmm

# Pruebo esto para que se use la misma versión de Perl y SignalP funcione 
PATH=/home/aredondo/miniconda3/envs/SqueezeMeta/bin:$PATH

# Prueba en un genoma
nohup python run_dbcan.py /media/disk2/aredondo/artico022_DNA/genomes/1347342.6/results/03.1347342.6.faa \
	protein --out_dir ./genome_output --hmm_eval 1e-15 --hmm_cpu 8 &

# Esto es el modo protein
nohup python run_dbcan.py /media/disk2/aredondo/artico022_DNA/results/03.artico022_DNA.faa \
        protein --out_dir ./output -c /media/disk2/aredondo/artico022_DNA/intermediate/03.artico022_DNA.gff \
        --hmm_eval 1e-15 --hmm_cpu 8 &

# Vamos a probarlo en modo meta
nohup python run_dbcan.py /media/disk2/aredondo/artico022_DNA/results/03.artico022_DNA.fna \
        meta --out_dir ./meta_output -c /media/disk2/aredondo/artico022_DNA/intermediate/03.artico022_DNA.gff \
        --hmm_eval 1e-15 --hmm_cpu 8 &

# Vamos a probarlo con los parámetros de SQM
nohup python run_dbcan.py /media/disk2/aredondo/artico022_DNA/results/03.artico022_DNA.fna \
        protein --out_dir ./sqm_par_output -c /media/disk2/aredondo/artico022_DNA/intermediate/03.artico022_DNA.gff \
        --dia_eval 1e-03 --hmm_eval 1e-10 --hmm_cpu 8 &

'''
        [inputFile] - FASTA format file of either nucleotide or protein sequences
                
        [inputType] - protein=proteome, prok=prokaryote, meta=metagenome/mRNA/CDSs/short DNA seqs
                
        [--out_dir] - REQUIRED, user specifies an output directory.
                
        [-c AuxillaryFile]- optional, include to enable CGCFinder. If using a proteome input,
        the AuxillaryFile must be a GFF or BED format file containing gene positioning
        information. Otherwise, the AuxillaryFile may be left blank.
                                        
        [-t Tools] - optional, allows user to select a combination of tools to run. The options are any
                                        combination of 'diamond', 'hmmer', and 'hotpep'. The default value is 'all' which runs all three tools.

	[--dia_eval] - optional, allows user to set the DIAMOND E Value. Default = 1e-121.
                
        [--dia_cpu] - optional, allows user to set how many CPU cores DIAMOND can use. Default = 5.
                
        [--hmm_eval] - optional, allows user to set the HMMER E Value. Default = 1e-35.
                
        [--hmm_cov] - optional, allows user to set the HMMER Coverage value. Default = 0.35.
                
        [--hmm_cpu] - optional, allows user to set how many CPU cores HMMER can use. Default = 1.
                
        [--hot_hits] - optional, allows user to set the Hotpep Hits value. Default = 4.
                
        [--hot_freq] - optional, allows user to set the Hotpep Frequency value. Default = 2.0.
                
        [--hot_cpu] - optional, allows user to set how many CPU cores Hotpep can use. Default = 4.
                
        [--out_pre] - optional, allows user to set a prefix for all output files.

		[--db_dir] - optional, allows user to specify a database directory. Default = db/
                
        [--cgc_dis] - optional, allows user to specify CGCFinder Distance value. Allowed values are integers between 0-10. Default = 2.
                
        [--cgc_sig_genes] - optional, allows user to specify CGCFinder Signature Genes. The options are, 'tp': TP and CAZymes, 'tf': TF and CAZymes, and 'all': TF, TP, and CAZymes. Default = 'tp'.

'''