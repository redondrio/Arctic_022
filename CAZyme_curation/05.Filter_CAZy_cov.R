#####################################
# FILTER CAZYMES BY CONSENSUS
# AND RECALCULATE COVERAGE
#####################################

# The idea here is to extract good CAZymes from the Arctic
# using the consensus criteria, and rescale those to
# trcn so that they can be analysed with the rest

# The idea is:
# 1. Import the reference annotations
# 2. Select cazymes with consensus
# 3. Aggregate by cazy family and sum coverage
# 4. Substitute in Arctic$functions$CAZy$cov

# The already saved object is here
# saveRDS(cazy_agg,file="filt_cazy_cov.rds")

cazy_agg <- readRDS("filt_cazy_cov.rds")

library(textclean)

fun2grp<-function(fun.name){
	fun.name<-gsub("\\*","",fun.name)
	grp.name<-fun.name
	lettersites<-grepl("^[A-Z]",fun.name)
	grp.name[lettersites]<-gsub("[_[0-9]]*","",fun.name)[lettersites]
	grp.name[!lettersites]<-gsub(
		"\\.[0-9]*\\.[0-9]*\\.[-[0-9]*]*$","",fun.name)[!lettersites]
	grp.name[grp.name==1]="AA"
	grp.name[grp.name==2]="GT"
	grp.name[grp.name==3]="EC3"
	grp.name[grp.name==4]="PL"
	return(grp.name)	
}

# 1. Import the reference annotations
setwd("/media/disk2/aredondo/artico022_DNA/genomes/CAZyDB/results/uniques")

gh_kegg <- read.table("CAZyDB.GH.curid.kegg")
gh_cog  <- read.table("CAZyDB.GH.curid.cog")
gh_pfam <- read.table("CAZyDB.GH.curid.pfam")

gt_kegg <- read.table("CAZyDB.GT.curid.kegg")
gt_cog  <- read.table("CAZyDB.GT.curid.cog")
gt_pfam <- read.table("CAZyDB.GT.curid.pfam")

pl_kegg <- read.table("CAZyDB.PL.curid.kegg")
pl_cog  <- read.table("CAZyDB.PL.curid.cog")
pl_pfam <- read.table("CAZyDB.PL.curid.pfam")

ce_kegg <- read.table("CAZyDB.CE.curid.kegg")
ce_cog  <- read.table("CAZyDB.CE.curid.cog")
ce_pfam <- read.table("CAZyDB.CE.curid.pfam")

cbm_kegg <- read.table("CAZyDB.CBM.curid.kegg")
cbm_cog  <- read.table("CAZyDB.CBM.curid.cog")
cbm_pfam <- read.table("CAZyDB.CBM.curid.pfam")

k_cazy <- list(gh_kegg,gt_kegg,pl_kegg,ce_kegg,cbm_kegg)
c_cazy <- list(gh_cog,gt_cog,pl_cog,ce_cog,cbm_cog)
p_cazy <- list(gh_pfam,gt_pfam,pl_pfam,ce_pfam,cbm_pfam)

filtered <- c("GH","GT","PL","CE","CBM")
names(k_cazy) <- filtered
names(c_cazy) <- filtered
names(p_cazy) <- filtered

# 2. Select cazymes with consensus
cons_vec <- readRDS("good_cazy_vec.rds")
cons_vec <- c()
orf_tab <- Arctic$orfs$table[Arctic$orfs$table$CAZy!="",c(9,12,15,16)]

for (orf_i in seq(1,nrow(orf_tab))){
	# Get the type of cazyme
	cons = 1
	if (orf_i%%10000 == 0){print(nrow(orf_tab)-orf_i)}
	cazy  <- fun2grp(orf_tab$CAZy[orf_i])
	#print(cazy)
	#print(orf_tab[orf_i,])
	if (cazy %in% filtered){
		# Get the annotations for that orf
		k_st <- orf_tab$"KEGG ID"[orf_i]
		c_st <- orf_tab$"COG ID"[orf_i]
		p_st <- orf_tab$PFAM[orf_i]
		k_id <- unlist(strsplit(k_st,"[;*]",))
		c_id <- unlist(strsplit(c_st,"[;*]",))
		p_id <- unlist(strsplit(p_st,"[; ]",))
		# Check if they are in the reference
		k_check <- any(k_id %in% unlist(k_cazy[cazy]))
		c_check <- any(c_id %in% unlist(c_cazy[cazy]))
		p_check <- any(p_id %in% unlist(p_cazy[cazy]))
		cons <- c(k_check,c_check,p_check)
		#print(cons)
	}
	# Add the tag
	if (sum(cons)>0 & cazy!=""){
		cons_vec <- c(cons_vec,TRUE)
	} else {
		cons_vec <- c(cons_vec,FALSE)
	}
}
# Save the object for future use
saveRDS(cons_vec,file="good_cazy_vec.rds")
cons_vec <- readRDS("good_cazy_vec.rds")

# Subset the coverage of orfs that are good cazymes
orf_ids <- rownames(Arctic$orfs$table[Arctic$orfs$table$CAZy!="",])[cons_vec]
cov_tab <- Arctic$orfs$table[orf_ids,c(16,42:65)]
cov_tab$CAZy <- gsub("\\*","",cov_tab$CAZy)

# Add excluded to check counts are conserved
exc_ids <- rownames(orf_tab)[!(cons_vec)]
unclass <- colSums(Arctic$orfs$table[exc_ids,c(42:65)])

# 3. Aggregate by cazy family and sum coverage
cazy_agg <- aggregate(cov_tab[,-1],list(cov_tab$CAZy),sum)
rownames(cazy_agg) <- cazy_agg[,1]
cazy_agg <- cazy_agg[,-1]
colnames(cazy_agg) <- colnames(Arctic$functions$CAZy$cov)

# Add unclassified reads
cazy_agg["Unclassified",] <- unclass+Arctic$functions$CAZy$cov["Unclassified",]

# 4. Substitute in Arctic$functions$CAZy$cov
Arctic$functions$CAZy$cov <- cazy_agg

# And substitute the CAZy annotation table with blank
Arctic$orfs$table$CAZy[!(rownames(
	Arctic$orfs$table) %in% orf_ids)] <- ""

# 5. Check if numbers make sense

# Check total sum is the same in the aggregated and original
# The aggregated has now moved many functions to Unclassified
colSums(cazy_agg)
colSums(Arctic$functions$CAZy$cov)

# Let's now check how much classified coverage was removed
# and how many orfs were unclassified
# Almost 60% (58.6%) of the orfs were removed, which on average
# contributed 63.1% of the coverage
# This is more or less what was removed in the genomes,

tot <- length(cons_vec) #total number of original cazyme orfs
fil <- sum(cons_vec) #total number of filtered cazyme orfs
print(1-fil/tot)

fil_cov <- colSums(cazy_agg[-nrow(cazy_agg),])
tot_cov <- colSums(Arctic$functions$CAZy$cov[-(
	nrow(Arctic$functions$CAZy$cov)),])
print(mean(1-fil_cov/tot_cov))

# Check both tables have the same number of rows
nrow(cazy_agg)
nrow(Arctic$functions$CAZy$cov)

# 72 families got lost due to the filter
# let's see which
agg_rows<-rownames(cazy_agg)
fun_rows<-rownames(Arctic$functions$CAZy$cov)
setdiff(fun_rows,agg_rows)

# Manual sum
colSums(Arctic$orfs$table[Arctic$orfs$table$CAZy=="GH23"
	| Arctic$orfs$table$CAZy=="GH23*",c(42:65)])
# Aggregated sum
cazy_agg["GH23",]
# Original
Arctic$functions$CAZy$cov["GH23",]






