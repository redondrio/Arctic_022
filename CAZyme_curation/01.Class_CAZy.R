####################################
# CLASSIFY CAZYMES BASED ON 
# CONSENSUS WITH REF & KCP
####################################

# This second part will try to generate a DF with columns
# This is what was done before having the CAZyDB annotations
# CAZY, COG, KEGG and REF to try and generate a classifier
# that can tell apart between properly annotated CAZymes and
# missanotations, based on the consensus with COG and KEGG

library(ggplot2)
library(reshape2)

genomeid_list <- c(
"1347342.6",
"376686.10",
"313598.6",
"1798225.3",
"1336794.4",
"1336795.4",
"1336804.3",
"2058137.3")

# This third part will try to generate a DF with columns
# This is what was done AFTER having the CAZyDB annotations
# CAZY, COG, KEGG and REF to try and generate a classifier
# that can tell apart between properly annotated CAZymes and
# missanotations, based on the consensus with CAZyDB annots
# This is the loop to be run in a more compact way so that it can be copy-pasted

# Function to read files safely, even when they're empty or do not exist
safe_read <- function(file_name){
	list_annot <- tryCatch(read.table(file=paste0(file_name),sep="\t",fill=T,quote=""),
		error = function(e) data.frame())
	return(list_annot)
}


once=TRUE
cazyfam="CE"

for (genomeid in genomeid_list){
	setwd(paste0("/media/disk2/aredondo/artico022_DNA/genomes/",genomeid,"/results/",cazyfam,"_quality"))
	# Load the IDs of reference and SQM-annotated CAZymes
	ref.ids <- tryCatch(read.table(file=paste0("ref_",cazyfam,"_ids.txt"),colClasses="character"),
		error = function(e) c())
	sqm.ids <- read.table(file=paste0("sqm_",cazyfam,"_ids.txt"),colClasses="character")
	all.ids <- unique(c(ref.ids[,1],sqm.ids[,1]))

	# Load the reference annotations obtained by annotating CAZy
	CAZyDB.dir  <-  "/media/disk2/aredondo/artico022_DNA/genomes/CAZyDB/results/uniques/" 
	consensus.cog  <- read.table(file=paste0(CAZyDB.dir,"CAZyDB.",cazyfam,".curid.cog"),colClasses="character",quote="",fill=T,header=FALSE)
	consensus.kegg <- read.table(file=paste0(CAZyDB.dir,"CAZyDB.",cazyfam,".curid.kegg"),colClasses="character",quote="",fill=T,header=FALSE)
	consensus.pfam <- read.table(file=paste0(CAZyDB.dir,"CAZyDB.",cazyfam,".curid.pfam"),colClasses="character",quote="",fill=T,header=FALSE)

	# Extract the related annotations
	TP.cog <- safe_read(paste0("TP_",cazyfam,"_cog.txt"))
	FP.cog <- safe_read(paste0("FP_",cazyfam,"_cog.txt"))
	FN.cog <- safe_read(paste0("FN_",cazyfam,"_cog.txt"))

	TP.kegg <- safe_read(paste0("TP_",cazyfam,"_kegg.txt"))
	FP.kegg <- safe_read(paste0("FP_",cazyfam,"_kegg.txt"))
	FN.kegg <- safe_read(paste0("FN_",cazyfam,"_kegg.txt"))
	
	TP.pfam <- safe_read(paste0("TP_",cazyfam,"_pfam.txt"))
	FP.pfam <- safe_read(paste0("FP_",cazyfam,"_pfam.txt"))
	FN.pfam <- safe_read(paste0("FN_",cazyfam,"_pfam.txt"))
	
	cog.names  <- rbind(TP.cog,FP.cog,FN.cog)
	kegg.names <- rbind(TP.kegg,FP.kegg,FN.kegg)
	pfam.names <- rbind(TP.pfam,FP.pfam,FN.pfam)

	# Get the IDs of the proteins that have those annotations
	cog.cons  <- cog.names[,1][cog.names[,2] %in% consensus.cog[,1]]
	kegg.cons <- kegg.names[,1][kegg.names[,2] %in% consensus.kegg[,1]]
	pfam.cons <- pfam.names[,1][pfam.names[,2] %in% consensus.pfam[,1]]

	# Create boolean vectors indicating if ID in lists
	label.ref <- sapply(all.ids,function(x) x %in% ref.ids[,1])
	label.sqm <- sapply(all.ids,function(x) x %in% sqm.ids[,1])
	label.cog <- sapply(all.ids,function(x) x %in% cog.cons)
	label.kegg <- sapply(all.ids,function(x) x %in% kegg.cons)
	label.pfam <- sapply(all.ids,function(x) x %in% pfam.cons)
	CAZyme.classifier <- data.frame(
		Ref=label.ref,
		SQM=label.sqm,
		COG=label.cog,
		KEGG=label.kegg,
		PFAM=label.pfam)

	# Add this so that all tables have the same shape
	CAZyme.classifier <- rbind(CAZyme.classifier,c(FALSE,FALSE,FALSE,FALSE,FALSE))
	CAZyme.classifier <- rbind(CAZyme.classifier,c(TRUE,TRUE,TRUE,TRUE,TRUE))
	# Print those that have a specific profile
	print(CAZyme.classifier[
		CAZyme.classifier$Ref==FALSE &
		CAZyme.classifier$SQM==TRUE &
		CAZyme.classifier$KEGG==FALSE &
		CAZyme.classifier$COG==TRUE | 
		CAZyme.classifier$PFAM==TRUE
		,])
	if (once){
		tab.clas <- table(CAZyme.classifier)
		once=FALSE
	} else {
		tab.clas <- tab.clas+table(CAZyme.classifier)
	}
	# Remove the all-positive that was added to guarantee the shape
	tab.clas[length(tab.clas)] <- tab.clas[length(tab.clas)]-1
}
print(tab.clas)






