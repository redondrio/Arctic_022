####################################
# EN R
# Este script repite el análisis que
# se ha hecho con el ártico pero con
# los genomas individuales
####################################

# Load filtered coverage table into Arctic
setwd("/media/disk2/aredondo/artico022_DNA/rds")
Arctic$functions$CAZy$cov<-readRDS("filt_cazy_cov.rds")

# Filtering starts gere
setwd(paste0("/media/disk2/aredondo/artico022_DNA/genomes/"))
genome_list<-read.table(".done_genomes.txt",sep="\t",colClasses="character")

# Manually curated associated domains 
GH.consensus.cogs<-read.table("./CAZyDB/results/uniques/CAZyDB.GH.curid.cog",header=FALSE)[,1]
GH.consensus.kegg<-read.table("./CAZyDB/results/uniques/CAZyDB.GH.curid.kegg",header=FALSE)[,1]
GH.consensus.pfam<-read.table("./CAZyDB/results/uniques/CAZyDB.GH.curid.pfam",header=FALSE)[,1]
GT.consensus.cogs<-read.table("./CAZyDB/results/uniques/CAZyDB.GT.curid.cog",header=FALSE)[,1]
GT.consensus.kegg<-read.table("./CAZyDB/results/uniques/CAZyDB.GT.curid.kegg",header=FALSE)[,1]
GT.consensus.pfam<-read.table("./CAZyDB/results/uniques/CAZyDB.GT.curid.pfam",header=FALSE)[,1]
PL.consensus.cogs<-read.table("./CAZyDB/results/uniques/CAZyDB.PL.curid.cog",header=FALSE)[,1]
PL.consensus.kegg<-read.table("./CAZyDB/results/uniques/CAZyDB.PL.curid.kegg",header=FALSE)[,1]
PL.consensus.pfam<-read.table("./CAZyDB/results/uniques/CAZyDB.PL.curid.pfam",header=FALSE)[,1]
CE.consensus.cogs<-read.table("./CAZyDB/results/uniques/CAZyDB.CE.curid.cog",header=FALSE)[,1]
CE.consensus.kegg<-read.table("./CAZyDB/results/uniques/CAZyDB.CE.curid.kegg",header=FALSE)[,1]
CE.consensus.pfam<-read.table("./CAZyDB/results/uniques/CAZyDB.CE.curid.pfam",header=FALSE)[,1]
CBM.consensus.cogs<-read.table("./CAZyDB/results/uniques/CAZyDB.CBM.curid.cog",header=FALSE)[,1]
CBM.consensus.kegg<-read.table("./CAZyDB/results/uniques/CAZyDB.CBM.curid.kegg",header=FALSE)[,1]
CBM.consensus.pfam<-read.table("./CAZyDB/results/uniques/CAZyDB.CBM.curid.pfam",header=FALSE)[,1]


genome.analysis<-data.frame(Avg.USi=c(),Pept=c(),Sulf=c(),
	SusC=c(),SusD=c(),Tandem=c(),Pept.to.CAZy=c(),CAZy=c(),GH=c(),
	GT=c(),PL=c(),CE=c(),AA=c(),CBM=c())
profile.analysis<-list()

find.tandem<-function(frame){
  #the input frame is expected to have two columns, index and homolog (C/D)
  #indices must be ordered
  tandem.vec<-c()
  for (i in 2:length(frame[,1])){
    if ((frame[i,1]-frame[(i-1),1]==1) & frame[i,2]!=frame[(i-1),2]){
      tandem.vec<-c(tandem.vec,frame[i-1,1],frame[i,1])
    } #the first condition was changed to ==2 to check for interrupted tandems (none were found)
  }
  return(tandem.vec)
}

fun2fam<-function(fun.name){
	fam.name<-gsub("\\.[-[0-9]*]*$","",fun.name)
	fam.name<-gsub("_[0-9]*$","",fam.name)
	return(fam.name)
}

for (genomecode in genome_list[,1]){
	# Set dir
	setwd(paste0("/media/disk2/aredondo/artico022_DNA/genomes/",genomecode,"/results/"))
	# Read tables
	print(genomecode)
	keggtable<-read.table(paste0("07.",genomecode,".fun3.kegg"),header=FALSE,fill=TRUE)
	cogtable<-read.table(paste0("07.",genomecode,".fun3.cog"),header=FALSE,fill=TRUE,quote="")
	pfamtable<-read.table(paste0("07.",genomecode,".fun3.pfam"),header=FALSE,sep="\t",quote="")
	cazytable<-read.table(paste0("07.",genomecode,".fun3.CAZyDB.dmnd"),header=FALSE,sep="\t",quote="")
	names(keggtable)<-c("ORF","BESTHIT","BESTAVER")
	names(cogtable)<-c("ORF","BESTHIT","BESTAVER")
	names(pfamtable)<-c("ORF","BESTHIT")
	names(cazytable)<-c("ORF","BESTHIT","BESTAVER")
	# Count USiCGs to check the assumption of single-copy
	usicgs.ids<-"PF00687|PF03946|PF00572|PF00238|PF00828|PF03947|PF00673|PF00347|PF00338|PF00411|PF00164|PF00416|PF00163|PF00177|PF00410"
	usicgs.list<-strsplit(usicgs.ids,"\\|")
	count=0
	for (usi in usicgs.list[[1]]){
		count<-count+length(grep(usi,pfamtable$BESTHIT))
	}
	avg.usi<-count/length(unlist(usicgs.list))
	# Count peptidase annotations
	pept.allindx<-as.numeric(grep("[P|p]eptidase|[P|p]rotease|[P|p]roteinase", pfamtable$BESTHIT))
	pept.escindx<-as.numeric(grep("[P|p]utative|[I|i]nhibitor", pfamtable$BESTHIT))
	pept.pepindx<-setdiff(pept.allindx,pept.escindx)
	pept.count<-length(pept.pepindx)
	# Count sulfatase annotations
	sulf.indx<-grep("Sulfatase",pfamtable$BESTHIT)
	sulf.count<-nrow(pfamtable[sulf.indx,])
	# Count Sus-homolog count
	susc.indx<-grep("PF00593|PF07715|PF13715",pfamtable$BESTHIT)
	susd.indx<-grep("PF07980|PF12741|PF12771|PF14322",pfamtable$BESTHIT)
	susc.count<-nrow(pfamtable[susc.indx,])
	susd.count<-nrow(pfamtable[susd.indx,])
	# Count SusCD tandems
	susc.frame<-data.frame(index=susc.indx,homolog=c("C"))
	susd.frame<-data.frame(index=susd.indx,homolog=c("D"))
	suscd.frame<-rbind(susc.frame,susd.frame)
	suscd.frame<-suscd.frame[order(suscd.frame$index),]
	tandems<-find.tandem(suscd.frame)
	tandem.count<-length(tandems)/2

	# CAZy profile
	cazy.count<-nrow(cazytable)
	GH.count<-nrow(cazytable[grep("GH",cazytable$BESTHIT),])
	GT.count<-nrow(cazytable[grep("GT",cazytable$BESTHIT),])
	PL.count<-nrow(cazytable[grep("PL",cazytable$BESTHIT),])
	CE.count<-nrow(cazytable[grep("CE",cazytable$BESTHIT),])
	AA.count<-nrow(cazytable[grep("AA",cazytable$BESTHIT),])
	CBM.count<-nrow(cazytable[grep("CBM",cazytable$BESTHIT),])

	filter <- function(fam,cons.k,cons.c,cons.p){
		# Filtered GH with COG, KEGG or PFAM consensus
		id.sqm <- cazytable[grep(fam,cazytable$BESTHIT),1]
		id.cons.cog <- cogtable[cogtable$BESTHIT %in% cons.c,1]
		id.cons.kegg <- keggtable[keggtable$BESTHIT %in% cons.k,1]
		pfam.vec <- levels(pfamtable$BESTHIT)[as.integer(pfamtable$BESTHIT)]
		pfam.list <- sapply(pfam.vec, function(x) strsplit(x," "))
		pfam.bool <- sapply(pfam.list, function(x) any(x %in% cons.p))
		id.cons.pfam <- pfamtable[pfam.bool,1]
		id.cons.cog <- levels(id.cons.cog)[as.integer(id.cons.cog)]
		id.cons.kegg <- levels(id.cons.kegg)[as.integer(id.cons.kegg)]
		id.cons.pfam <- levels(id.cons.pfam)[as.integer(id.cons.pfam)]
		cons.count <- length(id.sqm[id.sqm %in% c(id.cons.cog,id.cons.kegg,id.cons.pfam)])
		return(cons.count)
	}

	cons.GH.count <- filter("GH",GH.consensus.kegg,GH.consensus.cogs,GH.consensus.pfam)
	cons.GT.count <- filter("GT",GH.consensus.kegg,GH.consensus.cogs,GT.consensus.pfam)
	cons.PL.count <- filter("PL",PL.consensus.kegg,PL.consensus.cogs,PL.consensus.pfam)
	cons.CE.count <- filter("CE",CE.consensus.kegg,CE.consensus.cogs,CE.consensus.pfam)
	cons.CBM.count <- filter("CBM",CBM.consensus.kegg,CBM.consensus.cogs,CBM.consensus.pfam)

	# CAZyme profile
	cazytable$FAM<-fun2fam(cazytable$BESTHIT)
	cazy.profile<-table(cazytable$FAM)
	
	# Save all in a data frame and list
	genome.stats<-data.frame(avg.usi,pept.count,sulf.count,susc.count,
		susd.count,tandem.count,pept.count/cazy.count,cazy.count,
		cons.GH.count,cons.GT.count,cons.PL.count,cons.CE.count,
		cons.CBM.count,GT.count,PL.count,CE.count,AA.count,CBM.count)
	genome.analysis<-rbind(genome.analysis,genome.stats)
	profile.analysis[[genomecode]]<-cazy.profile
}
rownames(genome.analysis) <- genome_list[,1]
