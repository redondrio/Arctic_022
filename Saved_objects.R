##################################
#	SAVED FUNCTIONS
##################################

# Functions to clean CAZy names
# Just removes *
fun2clean <- function(fun.name){
	fun.name<-gsub("\\*","",fun.name)
	return(fun.name)
}

# Removes * and subfamilies
fun2fam <- function(fun.name){
	fun.name<-gsub("\\*","",fun.name)
	fam.name <- gsub("\\.[-[0-9]*]*$","",fun.name)
	fam.name <- gsub("_[0-9]*$","",fam.name)
	return(fam.name)
}

# Generates the big 6 groups
fun2grp <- function(fun.name){
	fun.name<-gsub("\\*","",fun.name)
	grp.name <- fun.name
	lettersites <- grepl("^[A-Z]",fun.name)
	grp.name[lettersites] <- gsub("[_[0-9]]*","",fun.name)[lettersites]
	grp.name[!lettersites] <- gsub(
		"\\.[0-9]*\\.[0-9]*\\.[-[0-9]*]*$","",fun.name)[!lettersites]
	grp.name[grp.name==1]="AA"
	grp.name[grp.name==2]="GT"
	grp.name[grp.name==3]="EC3"
	grp.name[grp.name==4]="PL"
	return(grp.name)	
}

# Generates the big 6 groups
fam2grp <- function(fam.name){
	fam.name<-gsub("\\*","",fam.name)
	grp.name <- fam.name
	lettersites <- grepl("^[A-Z]",fam.name)
	grp.name[lettersites] <- gsub("[0-9]*$","",fam.name)[lettersites]
	grp.name[!lettersites] <- gsub(
		"\\.[0-9]*\\.[0-9]*$","",fam.name)[!lettersites]
	grp.name[grp.name==1]="AA"
	grp.name[grp.name==2]="GT"
	grp.name[grp.name==3]="EC3"
	grp.name[grp.name==4]="PL"
	return(grp.name)	
}

# Function to generate the subsets with taxon-rescaled copy numbers
save.litesub<-function(bigsubset,name,rank,taxon,save=TRUE){
	print(paste0("Subsetting ",taxon))
	big_subset<-subsetTax(bigsubset,rank,taxon,rescale_tpm=T,rescale_copy_number=F)
	print("Reducing subset to lite")
	lite_subset<-combineSQMlite(big_subset)
	rm(big_subset)
	#Normalisation
	print("Normalising")
	eUSiPFAM <- USiPFAM[USiPFAM %in% rownames(lite_subset$functions$PFAM$cov)]
	USiPFAM_cov <- apply(lite_subset$functions$PFAM$cov[eUSiPFAM,], 2, median)
	PFAM_copynumber <- t(t(lite_subset$functions$PFAM$cov) / USiPFAM_cov)
	lite_subset$functions$PFAM$copy_number <- PFAM_copynumber

	eUSiCGs <- USiKEGG[USiKEGG %in% rownames(lite_subset$functions$KEGG$cov)]
	USiKEGG_cov <- apply(lite_subset$functions$KEGG$cov[eUSiCGs,], 2, median)
	KEGG_copynumber <- t(t(lite_subset$functions$KEGG$cov) / USiKEGG_cov)
	lite_subset$functions$KEGG$copy_number <- KEGG_copynumber
	# CAZymes are rescaled with KEGG coverages
	CAZy_copynumber <- t(t(lite_subset$functions$CAZy$cov) / USiKEGG_cov)
	lite_subset$functions$CAZy$copy_number <- CAZy_copynumber

	eUSiCOGs <- USiCOGs[USiCOGs %in% rownames(lite_subset$functions$COG$cov)]
	USiCOGs_cov <- apply(lite_subset$functions$COG$cov[eUSiCOGs,], 2, median)
	COG_copynumber <- t(t(lite_subset$functions$COG$cov) / USiCOGs_cov)
	lite_subset$functions$COG$copy_number <- COG_copynumber

	if(save){	
		print("Saving SQMlite")
		filename=paste0(name,"_sublite_norm.rds")
		saveRDS(lite_subset, file=filename)
		print(paste0("Done, SQMlite saved as: ",filename))
	}
	return(lite_subset)
}

# Function to generate de table for boxplots and time series
arrayTRCN<-function(sub_functions,sub_levels,
	aggr_fun,gen_fun="Functions",samples=DNAsamples,dates=DNAdoy,
	subsets=c(),subset_names=c()){
	if(length(sub_functions)!=length(sub_levels)){
		print("Different length in functions and levels")
		return()
	}
	if(length(sub_functions)!=length(aggr_fun)){
		print("Different length in functions and aggregates")
		return()
	}
	if(length(subsets)!=length(subset_names)){
		print("Different length in subsets and names")
		return()
	}
	trcn_array=c()
	for (subset_i in subsets){
		#trcn para las funciones de interes en el subset
		for (fun_i in seq(1,length(sub_functions))){
			fun_exist <- (sub_functions[[fun_i]] %in% rownames(subset_i[["functions"]][[sub_levels[fun_i]]][["copy_number"]]))		
			#print(fun_exist)
			if(sum(fun_exist)==0){
				trcn_array <- c(trcn_array,rep(0,length(samples)))
			}else{
				#print(sub_functions[[fun_i]][fun_exist])
				plotted<-plotFunctions(subset_i,fun_level=sub_levels[fun_i],
					samples=samples,count="copy_number",fun=sub_functions[[fun_i]][fun_exist])
				aggr<-aggregate(plotted$data$abun,list(plotted$data$sample),sum)
				#print(aggr)
				trcn_array <- c(trcn_array,aggr$x)
			}	
		}
	}
	dim(trcn_array)=c(length(samples),length(aggr_fun),length(subsets))
	dimnames(trcn_array)[[1]]=dates
	dimnames(trcn_array)[[2]]=aggr_fun
	dimnames(trcn_array)[[3]]=c(subset_names)
	trcn_array <- melt(trcn_array,id=dimnames(trcn_array)[[3]])
	colnames(trcn_array) <- c("Dates",gen_fun,"Taxon","trcn")
	return(trcn_array)
}

# Function to generate the dataframes for cazy profiles
aggr.cazyfam <- function(subset,name,samples=DNAsamples,quartile=0.6){
	# Subset taxon and samples
	cazy.profile <- as.data.frame(subset$functions$CAZy$copy_number[,c(samples)])
	# Get the mean
	trcn.mean <- rowMeans(cazy.profile)

	# May not want to aggregate subfamilies into families bc the SQM object doesn't,
	# each subfamily is counted independently in $functions
	f_agg <- aggregate(trcn.mean,
		list(fun2clean(names(trcn.mean))),sum)
	f_agg$taxon <- name

	# Prepare it for ggplot
	names(f_agg)=c("family","value","taxon")
	f_agg$group <- fun2grp(f_agg$family)
	f_agg <- subset(f_agg,family!="Unclassified")

	# Set quartile threshold and renames those
	f_agg$tag <- f_agg$family
	f_agg$tag[f_agg$value<quantile(f_agg$value[f_agg$group=="AA"],prob=quartile) & f_agg$group=="AA"] <- "Minor"
	f_agg$tag[f_agg$value<quantile(f_agg$value[f_agg$group=="GT"],prob=quartile) & f_agg$group=="GT"] <- "Minor"
	f_agg$tag[f_agg$value<quantile(f_agg$value[f_agg$group=="GH"],prob=quartile) & f_agg$group=="GH"] <- "Minor"
	f_agg$tag[f_agg$value<quantile(f_agg$value[f_agg$group=="PL"],prob=quartile) & f_agg$group=="PL"] <- "Minor"
	f_agg$tag[f_agg$value<quantile(f_agg$value[f_agg$group=="CE"],prob=quartile) & f_agg$group=="CE"] <- "Minor"
	f_agg$tag[f_agg$value<quantile(f_agg$value[f_agg$group=="EC3"],prob=quartile) & f_agg$group=="EC3"] <- "Minor"
	f_agg$tag[f_agg$value<quantile(f_agg$value[f_agg$group=="CBM"],prob=quartile) & f_agg$group=="CBM"] <- "Minor"

	return(f_agg)
}

save(fun2clean,fun2fam,fun2grp,fam2grp,
	save.litesub,arrayTRCN,aggr.cazyfam,
	file="utils.RData")

##################################
#	SAVED DOMAINS
##################################

# These are indices specific to the Arctic object for SusC and SusD
susc.allindx <- as.numeric(grep("PF07715|PF00593|PF13715|TIGR04056", Arctic$orfs$table$PFAM)) #not uptdated with TIGR04056
susd.allindx <- as.numeric(grep("PF07980|PF12741|PF12771|PF14322", Arctic$orfs$table$PFAM))
saveRDS(susc.allindx, file="susc_allindx.rds")
saveRDS(susd.allindx, file="susd_allindx.rds")

# These are the reference SusC/D domains for PFAM
ref.susc=c("PF00593 [TonB dependent receptor]",
            "PF07715 [TonB-dependent Receptor Plug Domain]",
            "PF13715 [CarboxypepD_reg-like domain]",
            "TIGR04056 []")
ref.susd=c("PF07980 [SusD family]",
           "PF14322 [Starch-binding associating with outer membrane]",
           "PF12741 [Susd and RagB outer membrane lipoprotein ]",
           "PF12771 [Starch-binding associating with outer membrane]")
other.references=c("PF00884 [Sulfatase]","PF14707 [C-terminal region of aryl-sulfatase]",
			"PF14292 [SusE outer membrane protein]",
			"PF16411 [Outer membrane protein SusF_SusE]")

# These are the reference SusC/D domains for KEGG
susc.KEGG=c("K03832","K21573")
susd.KEGG=c("K21572")

# These are the reference SusC/D domains for KEGG
susc.COG=c("COG0810") #actually it's for TonB
susd.COG=c("") # not found

# These are the USiCGs annotations in KEGG, COG and PFAM
# those that have 2 domains in PFAM are excluded because their
# counts are split in half between the domains, but probably
# not all have been annotated with both domains, so the coverage is useless
USiKEGG <- c(USiCGs)
USiCOGs <- c("COG0081","COG0080","COG0102","COG0093",
	"COG0200","COG0090","COG0094","COG0097",
	"COG0051","COG0100","COG0048","COG0099",
	"COG0522","COG0049","COG0096")
USiPFAM <- c("PF00687 [Ribosomal protein L1p/L10e family]",
	"PF03946 [Ribosomal protein L11, N-terminal domain]", # the RPL11 has two domains
	#"PF00298 [Ribosomal protein L11, RNA binding domain]",
	"PF00572 [Ribosomal protein L13]",
	"PF00238 [Ribosomal protein L14p/L23e]",
	"PF00828 [Ribosomal proteins 50S-L15, 50S-L18e, 60S-L27A]",
	"PF03947 [Ribosomal Proteins L2, C-terminal domain]", # the RPL2 has two dmoains
	#"PF00181 [Ribosomal Proteins L2, RNA binding domain]",
	"PF00673 [ribosomal L5P family C-terminus]", # the RPL5 has two domains
	#"PF00281 [Ribosomal protein L5]",
	"PF00347 [Ribosomal protein L6]",
	"PF00338 [Ribosomal protein S10p/S20e]",
	"PF00411 [Ribosomal protein S11]",
	"PF00164 [Ribosomal protein S12/S23]",
	"PF00416 [Ribosomal protein S13/S18]",
	"PF00163 [Ribosomal protein S4/S9 N-terminal domain]", # the RPS4 has two domains
	#"PF01479 [S4 domain]",
	"PF00177 [Ribosomal protein S7p/S5e]",
	"PF00410 [Ribosomal protein S8]")

# PEPTIDASE, SULFATASE AND CAZY DOMAINS

# For peptidases I did not search the IDs in the databases, but searched only those that were
# already in the Arctic. Peptidases are a wide family and selecting all IDs with peptidase in
# the name is useless.

# Peptidase PFAM domains
pept.allindx <- as.numeric(grep("[P|p]eptidase|[P|p]rotease|[P|p]roteinase", Arctic$orfs$table$PFAM))
pept.escindx <- as.numeric(grep("[P|p]utative|[I|i]nhibitor", Arctic$orfs$table$PFAM))
# This selects those that are not putative or inhibitors
pept.pepindx <- setdiff(pept.allindx,pept.escindx)
# This extracts all unique domains in those ORFs, be them peptidase or not
pept.domasoc <- unique(unlist(strsplit(unique(Arctic$orfs$table$PFAM[pept.pepindx]),";")))
# And here I select those that are peptidases again
pept.domindx <- grep("[P|p]eptidase|[P|p]rotease|[P|p]roteinase",pept.domasoc)
pept.domains <- sort(pept.domasoc[pept.domindx])

rm(pept.allindx, pept.domasoc, pept.escindx, pept.pepindx, pept.domindx)

# Peptidase KEGG domains
pept.allindx <- as.numeric(grep("[P|p]eptidase|[P|p]rotease|[P|p]roteinase", Arctic$orfs$table$KEGGFUN))
pept.escindx <- as.numeric(grep("[P|p]utative|[I|i]nhibitor", Arctic$orfs$table$KEGGFUN))
# This selects those that are not putative or inhibitors
pept.peptindx <- setdiff(pept.allindx,pept.escindx)
# This extracts all unique domains in those ORFs, be them peptidase or not
pept.domasoc <- unique(unlist(strsplit(unique(Arctic$orfs$table$"KEGG ID"[pept.peptindx]),";")))
# Here I remove the * for consensus annotation
pept.KEGG <- unlist(lapply(pept.domasoc,function(x) gsub("\\*","",x)))

rm(pept.allindx, pept.domasoc, pept.escindx, pept.peptindx, pept.domasoc)

# Peptidase COG domains
pept.allindx <- as.numeric(grep("[P|p]eptidase|[P|p]rotease|[P|p]roteinase", Arctic$orfs$table$COGFUN))
pept.escindx <- as.numeric(grep("[P|p]utative|[I|i]nhibitor", Arctic$orfs$table$COGFUN))
# This selects those that are not putative or inhibitors
pept.peptindx <- setdiff(pept.allindx,pept.escindx)
# This extracts all unique domains in those ORFs, be them peptidase or not
pept.domasoc <- unique(unlist(strsplit(unique(Arctic$orfs$table$"COG ID"[pept.peptindx]),";")))
# Here I remove the * for consensus annotation
pept.COG <- unlist(lapply(pept.domasoc,function(x) gsub("\\*","",x)))

rm(pept.allindx, pept.domasoc, pept.escindx, pept.peptindx, pept.domasoc)

# Sulfatase
# The PFAM sulfatase is included in other.references

# Extracted from KEGG searching sulfatase
sulf.KEGG <- c("K01130","K01131","K01132","K01133","K01134","K01135",
	"K01136","K01137","K01138","K12374","K12375","K12376","K12381",
	"K14607","K18222","K22295","K22966","K24118")
# Extracted from COG searching sulfatase
sulf.COG <- c("COG2015","COG3119")

# CAZymes
cazy.orfs <- rownames(Arctic$orfs$table[Arctic$orfs$table$CAZy!="",])
GH.orfs <- rownames(Arctic$orfs$table)[grep("GH",Arctic$orfs$table$CAZy)]
GT.orfs <- rownames(Arctic$orfs$table)[grep("GT",Arctic$orfs$table$CAZy)]
PL.orfs <- rownames(Arctic$orfs$table)[grep("PL",Arctic$orfs$table$CAZy)]
CE.orfs <- rownames(Arctic$orfs$table)[grep("CE",Arctic$orfs$table$CAZy)]
AA.orfs <- rownames(Arctic$orfs$table)[grep("AA",Arctic$orfs$table$CAZy)]
CBM.orfs <- rownames(Arctic$orfs$table)[grep("CBM",Arctic$orfs$table$CAZy)]

# CAZyme subgroups
# There are 6 groups at CAZy
# Glycoside Hydrolases 	 (GHs) : hydrolysis and/or rearrangement of glycosidic bonds (see CAZypedia definition)
# Glycosyl Transferases  (GTs) : formation of glycosidic bonds (see definition)
# Polysaccharide Lyases  (PLs) : non-hydrolytic cleavage of glycosidic bonds
# Carbohydrate Esterases (CEs) : hydrolysis of carbohydrate esters
# Auxiliary Activities	 (AAs) : redox enzymes that act in conjunction with CAZymes.
# Carbohydrate Binding 	 (CBM) : carbohydrate binding modules

cazy.domains <- unique(unlist(strsplit(unique(Arctic$orfs$table[cazy.orfs,"CAZy"]),";")))
GH.domains <- unique(cazy.domains[grep("GH",cazy.domains)])
GT.domains <- unique(cazy.domains[grep("GT",cazy.domains)])
PL.domains <- unique(cazy.domains[grep("PL",cazy.domains)])
CE.domains <- unique(cazy.domains[grep("CE",cazy.domains)])
AA.domains <- unique(cazy.domains[grep("AA",cazy.domains)])
CBM.domains <- unique(cazy.domains[grep("CBM",cazy.domains)])

# Aquí se guarda el objeto
save(ref.susc,ref.susd,other.references,susc.KEGG,susd.KEGG,susc.COG,susd.COG,
	USiCOGs,USiPFAM,USiKEGG,sulf.KEGG,sulf.COG,pept.KEGG,pept.COG,
	pept.domains,cazy.domains,GH.domains,GT.domains,PL.domains,
	CE.domains,AA.domains,CBM.domains,file="Domains.RData")

##################################
#	SAVED SAMPLE DATA
##################################

DNAsamples <- c("9_Mar","13_Mar","17_Mar","23_Apr","1_May","5_May","10_May",
			"19_May","1_Jun","11_Jun","15_Jun","23_Jun","30_Jul")
RNAsamples <- c("9_Mar_RNA","11_Mar_RNA","17_Mar_RNA","23_Apr_RNA","1_May_RNA","15_May_RNA",
			"19_May_RNA","1_Jun_RNA","11_Jun_RNA","15_Jun_RNA","23_Jun_RNA")
goodRNAsamples <- c("11_Mar_RNA","17_Mar_RNA","15_May_RNA",
			"1_Jun_RNA","11_Jun_RNA","15_Jun_RNA","23_Jun_RNA")

DNAdates <- as.Date(c("2000-03-09","2000-03-13","2000-03-17","2000-04-23","2000-05-01","2000-05-05",
	"2000-05-10","2000-05-19","2000-06-01","2000-06-11","2000-06-15","2000-06-23","2000-07-30"))
RNAdates <- as.Date(c("2000-03-09","2000-03-11","2000-03-17","2000-04-23","2000-05-01","2000-05-15",
	"2000-05-19","2000-06-01","2000-06-11","2000-06-15","2000-06-23"))
goodRNAdates <- as.Date(c("2000-03-11","2000-03-17","2000-05-15",
	"2000-06-01","2000-06-11","2000-06-15","2000-06-23"))

DNAdoy <- c(68,72,76,113,121,125,130,139,152,162,166,174,211)
RNAdoy <- c(68,70,76,113,121,135,139,152,162,166,174)
goodRNAdoy <- c(70,76,135,152,162,166,174)

taxons <- c("superkingdom","phylum","class","order","family","genus","species")
totalMb <- sum(Arctic$contigs$table$Length)

save(DNAsamples,RNAsamples,goodRNAsamples,
	DNAdates,RNAdates,goodRNAdates,DNAdoy,RNAdoy,goodRNAdoy,
	taxons,totalMb,
	file="sample_data.RData")



