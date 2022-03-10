####################################################################
####################################################################
#	SUSC-SUSD TANDEMS
####################################################################
####################################################################

## Selected indices for subsetting rows
## 11210 genes with SusC domains, 1415 with SusD
susc.allindx <- as.numeric(grep("PF00593|PF07715|PF13715|TIGR04056", Arctic$orfs$table$PFAM))
susd.allindx <- as.numeric(grep("PF07980|PF12741|PF12771|PF14322", Arctic$orfs$table$PFAM))

## Taxonomic association of SusC and SusD domains

## Out of 9725 tax-classified genes with SusC domains,
## 		5279 are classified as Bacteroidetes (54.3%)
##		4209 are classified as Proteobacteria (43.3%)
##			of which 2979 are Gammaproteobacteria (30.6%)
##
## The PF13715 domain is almost exclusively from Bacteroidetes 
## (96.8% of classified) so I excluded ORFs with only that domain
## because they are not enoguh to represent a SusC homolog
## Out of 7786 tax-classified genes with SusC domains, excluding PF13715-only,
## 		3404 are classified as Bacteroidetes (43.7%)
##		4200 are classified as Proteobacteria (53.9%)
##			of which 2976 are Gammaproteobacteria (38.2%)

susc.longtax <- Arctic$orfs$table$Tax[susc.allindx]
susc.phytax <- table(unlist(sapply(susc.longtax, function(x)
	unlist(strsplit(x,";"))[grep("p_",unlist(strsplit(x,";")))])))
susc.clatax <- table(unlist(sapply(susc.longtax, function(x)
	unlist(strsplit(x,";"))[grep("c_",unlist(strsplit(x,";")))])))
susc.phylist <- unlist(sapply(susc.longtax, function(x)
	unlist(strsplit(x,";"))[grep("p_",unlist(strsplit(x,";")))]))

## Out of 1415 genes with SusD domains,
## 1300 are classified as Bacteroidetes (91.9%)

susd.longtax <- Arctic$orfs$table$Tax[susd.allindx]
susd.phytax <- table(unlist(sapply(susd.longtax, function(x)
	unlist(strsplit(x,";"))[grep("p_",unlist(strsplit(x,";")))])))
rm(susc.longtax)
rm(susd.longtax)

## Generate the frame for the function
susc.frame <- data.frame(index=susc.allindx,homolog=c("C"))
susd.frame <- data.frame(index=susd.allindx,homolog=c("D"))
suscd.frame <- rbind(susc.frame,susd.frame)
suscd.frame <- suscd.frame[order(suscd.frame$index),]

## Function to extract contiguous susC-susD indices
## I  checked that there were no tandems with a gene in between (there were none)
## and that there were no split susC genes (there was 1)

find.tandem <- function(frame){
  #the input frame is expected to have two columns, index and homolog (C/D)
  #indices must be ordered
  tandem.vec <- c()
  for (i in 2:length(frame[,1])){
    if ((frame[i,1]-frame[(i-1),1]==1) & frame[i,2]!=frame[(i-1),2]){
      tandem.vec <- c(tandem.vec,frame[i-1,1],frame[i,1])
    } #the first condition was changed to ==2 to check for interrupted tandems (none were found)
  }
  return(tandem.vec)
}

tandem.allindx <- find.tandem(suscd.frame) #362 pairs

# Classify tandem candidates based on domains
## This function generates a dataframe with most of the needed information of the candidate tandems

clas.tandem <- function(tandem.indx,tandem.annot,ref.susc,ref.susd){
  tf <- data.frame(index1=c(),index2=c(),matchC=c(),matchD=c(),otherC=c(),otherD=c(),CD=c())
  ref.pfam=c(ref.susc,ref.susd)
  for(i in seq(from=1,to=length(tandem.indx),by=2)){
    #get the domains as lists
    can1.pfam <- strsplit(tandem.annot[i],split=";")[[1]]
    can2.pfam <- strsplit(tandem.annot[i+1],split=";")[[1]]
    #if ever wanted to include the actual domain name, change toString for can1.pfam[which(...)]
    match1.pfam <- toString(which(ref.pfam %in% can1.pfam)) #ref indices to matches
    other1.pfam <- toString(which(!(can1.pfam %in% ref.pfam))) #candidate indices to others
    match2.pfam <- toString(which(ref.pfam %in% can2.pfam))
    other2.pfam <- toString(which(!(can2.pfam %in% ref.pfam)))
    CDorder=any(can1.pfam %in% ref.susc) #True: CD, False: DC
    if (CDorder){
      candidate <- data.frame(indexC=tandem.indx[i],
                      indexD=tandem.indx[i+1],
                      matchC=match1.pfam,
                      matchD=match2.pfam,
                      otherC=other1.pfam,
                      otherD=other2.pfam,
                      CD=CDorder)
    }else{
      candidate <- data.frame(indexC=tandem.indx[i+1],
                      indexD=tandem.indx[i],
                      matchC=match2.pfam,
                      matchD=match1.pfam,
                      otherC=other2.pfam,
                      otherD=other1.pfam,
                      CD=CDorder)
    }
    #print(candidate)
    tf <- rbind(tf,candidate)
  }
  return(tf)
}

tf <- clas.tandem(tandem.allindx,Arctic$orfs$table$PFAM[tandem.allindx],ref.susc,ref.susd)
rm(susc.frame)
rm(susd.frame)
rm(tandem.allindx)

'''
OLDER VERSION
clas.tandem <- function(tandem.indx,tandem.annot,ref.pfam,ref.susc){
  tf <- data.frame(index1=c(),index2=c(),matchC=c(),matchD=c(),otherC=c(),otherD=c(),CD=c())
  for(i in seq(from=1,to=length(tandem.indx),by=2)){
    #get the domains as lists
    can1.pfam <- strsplit(tandem.annot[i],split=";")[[1]]
    can2.pfam <- strsplit(tandem.annot[i+1],split=";")[[1]]
    matchC.pfam <- toString(which(ref.pfam %in% can1.pfam)) #ref indices to matches
    otherC.pfam <- toString(which(!(can1.pfam %in% ref.pfam))) #candidate indices to others
    matchD.pfam <- toString(which(ref.pfam %in% can2.pfam))
    otherD.pfam <- toString(which(!(can2.pfam %in% ref.pfam)))
    CDorder=any(can1.pfam %in% ref.susc) #True: CD, False: DC, if false, swap values
    if (! CDorder){temp <- matchC.pfam; matchC.pfam <- matchD.pfam; matchD.pfam <- temp}
    if (! CDorder){temp <- otherC.pfam; otherC.pfam <- otherD.pfam; otherD.pfam <- temp}
    #print(tandem.annot[i]); print(tandem.annot[i+1])
    #readline(prompt="Continue:")
    candidate <- data.frame(index1=tandem.indx[i],
                          index2=tandem.indx[i+1],
                          matchC=matchC.pfam,
                          matchD=matchD.pfam,
                          otherC=otherC.pfam,
                          otherD=otherD.pfam,
                          CD=CDorder)
    #print(candidate)
    #readline(prompt="Next:")
    tf <- rbind(tf,candidate)
  }
  return(tf)
}
'''

## Add gene names and contigs
tf$namesC <- gsub("\\.1","",rownames(Arctic$orfs$table)[tf$indexC])
tf$namesD <- gsub("\\.1","",rownames(Arctic$orfs$table)[tf$indexD])
tf$contigC <- Arctic$orfs$table$"Contig ID"[tf$indexC]
tf$contigD <- Arctic$orfs$table$"Contig ID"[tf$indexD]

## Extract positions and calculate distance
tf$strC <- sapply(tf$namesC,function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][1]))
tf$endC <- sapply(tf$namesC,function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][2]))
tf$strD <- sapply(tf$namesD,function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][1]))
tf$endD <- sapply(tf$namesD,function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][2]))
tf$dist <- 0
tf$dist[tf$CD==TRUE]=tf$strD[tf$CD==TRUE]-tf$endC[tf$CD==TRUE]
tf$dist[tf$CD==FALSE]=tf$strC[tf$CD==FALSE]-tf$endD[tf$CD==FALSE]

## Taxonomic classification 
tf$contig_phy <- Arctic$contigs$tax[tf$contigC,2]
#tf$contig_cla <- Arctic$contigs$tax[tf$contigC,3]
#tf$contig_ord <- Arctic$contigs$tax[tf$contigC,4]
#tf$contig_fam <- Arctic$contigs$tax[tf$contigC,5]
#tf$contig_gen <- Arctic$contigs$tax[tf$contigC,6]
#tf$contig_spe <- Arctic$contigs$tax[tf$contigC,7]

# So far we have 362 tandems


####################################################################
####################################################################
#	CURATING TANDEMS
####################################################################
####################################################################

# Remove tandems that are not in the same contig (1 tandem removed)
tf <- subset(tf,contigC==contigD)
tf <- subset(tf,select=-contigD)
colnames(tf)[colnames(tf)=="contigC"] <- "contig"

## Add tandem names
rownames(tf) <- sprintf("tandem_%s",seq(1:nrow(tf)))

# Extract indices of tandems in the same contig
mult_tandem <- duplicated(tf$contig)|duplicated(tf$contig,fromLast=T)

## Check for clusters of susCD genes that have been marked as separate tandems
## This was done by hand by checking tf[mult_tandem,] and selecting only those tandems that overlapped
## e.g. 1 and 2 were two independent tandems on the same contig, so they were not excluded,
## but 12 and 13 were two overlapping tandems, so it was actually a single 3-gene cluster
## manual_mult_tandem <- c(12,13,84,85,86,219,220,232,233,234,246,247,317,318,319)
## After annotating with TIGRfam the indices changed
manual_mult_tandem <- c(4,5,15,16,97,98,99,112,113,259,260,274,275,276,
  291,292,333,334,343,344,367,368,369,385,386,387)

## After getting the clusters, I checked for the intergenic distance and removed the cluster where genes were farther apart.
## In 4-gene clusters, the middle tandem was removed, leaving two continuous tandems
## e.g. in the 12-13 tandems, genes in tandem 12 were 12pb apart, but 289pb apart in tandem 13, so tandem 13 was dismissed.
##  GENE1---GENE2-------------------------GENE3
## |--Tndm12-|-|------------Tndm13-------------|
manual_false_tandem <- c(4,16,98,112,260,275,
  292,333,343,368,386)
# Remove tandems that overlap (9 tandems removed)
tf <- tf[-manual_false_tandem,]

# Extract contigs that may have an undetected PF00593 domain
blasters <- tf[(tf$matchC=='2' | tf$matchC=='3' | tf$matchC=='2, 3'),"contig"]
blast.contig.seq <- Arctic$contigs$seqs[blasters]
write.table(cbind(names(blast.contig.seq),blast.contig.seq),"blasters.fa",sep="\n",col.names=F,row.names=F)
# awk '/"m/ { print ">"$0 } ; /"G|C|T|A/ { gsub(/"/, "") ; print }' blasters.fa > blasters.fasta

## Add tandem names again
rownames(tf) <- sprintf("tandem_%s",seq(1:nrow(tf)))



####################################################################
####################################################################
#	CAZYMES AND PFAM IN TANDEM
####################################################################
####################################################################

## Extract cazymes/PFAM in the same contigs as the tandems
other_PFAM <- Arctic$orfs$table$PFAM[Arctic$orfs$table$"Contig ID" %in% tf$contig]
names(other_PFAM) <- Arctic$orfs$table$"Contig ID"[Arctic$orfs$table$"Contig ID" %in% tf$contig]
other_CAZy <- Arctic$orfs$table$CAZy[Arctic$orfs$table$"Contig ID" %in% tf$contig]
names(other_CAZy) <- Arctic$orfs$table$"Contig ID"[Arctic$orfs$table$"Contig ID" %in% tf$contig]

## Subset to eliminate blanks, SusC/D domains and unknown function domains
other_PFAM <- subset(other_PFAM,other_PFAM!="")
other_PFAM <- subset(other_PFAM,sapply(other_PFAM,function(x) !any(
	grep("PF00593|PF07715|PF13715|PF07980|PF12741|PF12771|PF14322|unknown",x))))
other_CAZy <- subset(other_CAZy,other_CAZy!="")

## Sum the extra annotations
sum_PFAM <- aggregate(other_PFAM,list(names(other_PFAM)),length) #232 contigs with other domains
sum_CAZy <- aggregate(other_CAZy,list(names(other_CAZy)),length) #201 contigs with CAZymes
sum_both <- tapply(
	c(sum_PFAM$x,sum_CAZy$x),
	c(sum_PFAM$Group.1,sum_CAZy$Group.1),sum)

## Add the number of extra annotations to each tandem
tf$extra <- sum_both[match(tf$contig,names(sum_both))]
tf$extra[is.na(tf$extra)] <- 0

## Get surrounding genes
## Modified to use only the susD gene as a reference of the PUL
nearby.dist <- function(margin){
	tndmlist <- list()
	for (contig in unique(tf$contig)){
		for (tndm_name in rownames(tf)[tf$contig==contig]){
      # Get all orfs in the same contig
			pul.can <- rownames(Arctic$orfs$table)[Arctic$orfs$table$"Contig ID"==contig]
			# Get the position of the SusD gene in the contig
			pul.dindx <- which(pul.can==tf[tndm_name,"namesD"])
      # Extract the range of genes to check
      range <- c(
				pul.dindx-margin,
				pul.dindx+margin)
			if(range[1]<0){range[1]=1}
			if(range[2]>length(pul.can)){range[2]=length(pul.can)}
			range <- seq(range[1],range[2])
			pul.can <- pul.can[range]
      # Get the start and end positions of each gene
			starts <- sapply(pul.can,function(x) round(
        as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][1]),0))
			ends <- sapply(pul.can,function(x) round(
        as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][2]),0))
			starts <- c(starts)
			ends <- c(0,ends[-length(ends)])
			# Calculate distances
      distances <- starts-ends
      # Get PFAM annotations
      p_annots <- c()
      c_annots <- c()
      for (gene in pul.can){
        p_annots <- c(p_annots, Arctic$orfs$table[gene,"PFAM"])
        c_annots <- c(c_annots, Arctic$orfs$table[gene,"CAZy"])
      }
      tndm_add <- data.frame(distances,p_annots,c_annots)
      rownames(tndm_add) <- pul.can
      tndmlist[[tndm_name]] <- tndm_add
      print(tndm_name)
		}
	}
	return(tndmlist)
}

distances <- nearby.dist(10)
saveRDS(distances,file="distance10.rds")
distance10 <- readRDS("distance10.rds")

# Export to .csv to load in python
for(i in seq(1,length(distances))){
  df <- distances[[i]]
  print(i)
  write.table(df, file="PULs.csv", append=TRUE, sep = ",")
}
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/rds/PULs.csv /Users/alvaro/Desktop/CNB-Proyecto/Rproject/PULs.csv

##################################
# TANDEM FRAME SAVED HERE
##################################

saveRDS(tf,file="tandem_frame.rds") #429 tandems in 416 contigs

##################################
#	GENERAL VIEW
##################################

setwd("/media/disk2/aredondo/artico022_DNA/graphs")

# Checking for genes too far apart to be the same tandem (none were excessively apart)
pdf(file="distplot_cur.pdf")
hist <- hist(tf$dist,ylim=c(0,20))
text(hist$mids,rep(15,length(hist$mids)),labels=hist$counts,adj=c(0.5,-0.5))
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/distplot_cur.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/distplot_cur.pdf

# Classifying the tandems by domain architecture
tandem.archs <- ftable(xtabs(~ matchC + matchD, data=tf))

# Checking for possibly truncated genes
pdf(file='lenplot_cur.pdf')
plot(tf$strC,tf$endC-tf$strC,log = 'x',xlab="susC start position",ylab='susC length (pb)',main ='Length of susC candidates')
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/lenplot_cur.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/lenplot_cur.pdf

## The lenplot suggests that genes lacking the 2,3 domains are possibly truncated,
## as most of them are found at the beginning of the contig
aggregate(tf$strC[tf$CD==TRUE],list(tf$matchC[tf$CD==TRUE]),function(x) sum(x<10))

####################################################################
####################################################################
# CHECKING EXCLUDED SUSC/D
####################################################################
####################################################################

# The aim here is to have a look at the SusC and SusD that are not included
# in a tandem and have thus been excluded from the analysis
# The main hypothesis is that these genes will be at the ends of contigs
# and will likely be part of an incomplete tandem that was not fully retrieved

# This is the test for the main hypothesis
# The idea is to compare the start and end positions of the genes
# with the start and end of the contig and see how much is left

# Extract the orfs with SusC that are not in a tandem
susc.allnames   <-  rownames(Arctic$orfs$table)[susc.allindx]
susc.lonenames  <-  susc.allnames[!(susc.allnames %in% tf$namesC)]
# Get their start and end positions
susc.start  <-  sapply(susc.lonenames,
  function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][1]))
susc.end    <-  sapply(susc.lonenames,
  function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][2]))
# Get their contigs and their length



contigcL <- Arctic$contigs$table[susc.lonenames,"Length"]
# Calculate the distance to the end of the contig
susc.left <- contigcL-susc.end #bases left from the end of the SusC to the end of the contig

susc.contigs <- unique(Arctic$orfs$table$"Contig ID"[susc.allindx])
susc.contigs <- subset(susc.contigs,!(susc.contigs %in% tf$contig)) #contig ID
susc.again <- as.numeric(grep("PF00593|PF07715|PF13715",
  Arctic$orfs$table$PFAM[Arctic$orfs$table$"Contig ID" %in% susc.contigs])) #contig index
susc.start <- sapply(rownames(Arctic$orfs$table)[susc.again],
  function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][1]))
susc.end <- sapply(rownames(Arctic$orfs$table)[susc.again],
  function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][2]))
contigcL <- Arctic$contigs$table[Arctic$orfs$table$"Contig ID"[susc.again],"Length"]
susc.left <- contigcL-susc.end #bases left from the end of the SusC to the end of the contig

# Same for SusD, extract those not in a tandem
susd.orfs <- rownames(Arctic$orfs$table)[susd.allindx]
susd.contigs <- unique(Arctic$orfs$table[susd.orfs,"Contig ID"])
susd.lone.contigs <- subset(susd.contigs,!(susd.contigs %in% tf$contig)) #contig ID
susd.again <- as.numeric(grep("PF07980|PF12741|PF12771|PF14322", 
  Arctic$orfs$table$PFAM[Arctic$orfs$table$"Contig ID" %in% susd.lone.contigs])) #contig index
susd.start <- sapply(rownames(Arctic$orfs$table)[susd.again],
  function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][1]))
susd.end <- sapply(rownames(Arctic$orfs$table)[susd.again],
  function(x) as.numeric(strsplit(strsplit(x,'_')[[1]][3],'-')[[1]][2]))
contigdL <- Arctic$contigs$table[Arctic$orfs$table$"Contig ID"[susd.again],"Length"]
susd.left <- contigdL-susd.end #bases left from the end of the SusD to the end of the contig

# Prepare the data for plotting
plotting <- data.frame(rbind(
  cbind(
    cbind(susc.start,susc.left),
    "SusC"),
  cbind(
    cbind(susd.start,susd.left),
    "SusD"))
  )
colnames(plotting) <- c("Start","Left","Homolog")
plotting[,1] <- as.numeric(as.character(plotting[,1]))
plotting[,2] <- as.numeric(as.character(plotting[,2]))
plotting[plotting==0] <- 0.1 #to avoid problems with the log scale

# Plotting
# The plot includes lines to represent the avg length of SusC and SusD
# to see if the observed genes have enough "space" to accomodate a SusC/SusD gene
# The hypothesis is that they won't have enough space, supporting the idea that
# the tandem is incomplete
pdf(file="positions_sus.pdf")
ggplot(plotting, aes(x=Start,y=Left,group=Homolog)) +
  geom_point(aes(color=Homolog)) +
  scale_x_continuous(trans = 'log2') + #log scale
  scale_y_continuous(trans = 'log2') + #log scale
  geom_hline(yintercept=412, color="red") + geom_vline(xintercept=412, color="red") + #avg SusD length
  geom_hline(yintercept=1280, color="darkcyan") + geom_vline(xintercept=1280, color="darkcyan") #avg SusC length
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/positions_sus.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/positions_sus.pdf

middle <- plotting[plotting$Start>412,]
middle <- middle[middle$Homolog=="SusC",]
middle <- middle[middle$Left>412,]
print(middle)
middle <- plotting[plotting$Start>1280,]
middle <- middle[middle$Homolog=="SusD",]
middle <- middle[middle$Left>1280,]
print(middle)
# Out of 10837 SusC homologs, 1058 (10%) have space for a SusD at both sides
# Out of 1054 SusD homologs, 21 (2%) have space for a SusC at both sides

# Another hypothesis is that the SusC/D not in tandem have other functions
# The idea here is to see if these SusC/D are associated with genes involved in
# activies other than glycan metabolism, such as metal ion uptake, by looking
# at the functions of the associated domains in those contigs
susc.subset <- subsetContigs(Arctic,susc.contigs)
pdf(file="suscassoc.pdf")
PFAMfun <- plotFunctions(susc.subset,fun_level="PFAM",samples=DNAsamples)
print(PFAMfun)
KEGGfun <- plotFunctions(susc.subset,fun_level="KEGG",samples=DNAsamples)
print(KEGGfun)
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/suscassoc.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/suscassoc.pdf
susd.subset <- subsetContigs(Arctic,susd.contigs)
pdf(file="susdassoc.pdf")
PFAMfun <- plotFunctions(susd.subset,fun_level="PFAM",samples=DNAsamples)
print(PFAMfun)
KEGGfun <- plotFunctions(susd.subset,fun_level="KEGG",samples=DNAsamples)
print(KEGGfun)
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/susdassoc.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/susdassoc.pdf


# Extracting sequences that may be undetected SusC
susd.orfs <- rownames(Arctic$orfs$table)[susd.allindx]
susd.contigs <- unique(Arctic$orfs$table[susd.orfs,"Contig ID"])
susd.lone.contigs <- subset(susd.contigs,!(susd.contigs %in% tf$contig)) #contig ID
lone.assoc.orfs <- Arctic$orfs$table[Arctic$orfs$table$"Contig ID" %in% susd.lone.contigs,c("PFAM","Contig ID")]
blank.lone.orfs <- rownames(lone.assoc.orfs[lone.assoc.orfs$PFAM=="",])
hmmers <- Arctic$orfs$seqs[blank.lone.orfs]
write.table(hmmers,file="hmmers_undet_SusC.txt")

cat hmmers_undet_SusC.txt | tr "\"\ \"" "\n"| tr "mega" ">mega" | grep . > hmmers_undet_SusC.fa
cat hmmers_undet_SusC.fa | tr -d \* > hmmers_undet_SusC.fasta
# nano to remove the first x

# HMMer search
wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz
tar -xzf hmm_PGAP.HMM.tgz
cat hmm_PGAP/* >> tigrfam.hmm
hmmpress tigrfam.hmm
hmmscan --tblout hmm_table.out --noali tigrfam.hmm hmmers_undet_SusC.fasta




