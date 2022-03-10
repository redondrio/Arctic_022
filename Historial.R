Socorro at ngarcia@cnb.csic.es

##################################
#	CONNECTING
##################################
'''
ext-alvaro.redondo
ssh -y aredondo@zobel.cnb.csic.es
'''

## CONDA ACTIVATE SQM1_4
library("SQMtools")
setwd("/media/disk2/aredondo")
## CONDA ACTIVATE SQM1_4
Arctic <- loadSQM("artico022_DNA", engine="data.table")

setwd("/media/disk2/aredondo/artico022_DNA/rds")
load("Domains.RData")
load("sample_data.RData")
load("utils.RData")
library(ggplot2)
library(reshape)
library(gdata)
tf <- readRDS("tandem_frame.rds")
susc.allindx <- readRDS("susc_allindx.rds")
susd.allindx <- readRDS("susd_allindx.rds")
rdsdir <- "/media/disk2/aredondo/artico022_DNA/rds"
graphdir <- "/media/disk2/aredondo/artico022_DNA/graphs"

cazy_agg <- readRDS("filt_cazy_cov.rds")
Arctic$functions$CAZy$cov <- as.matrix(cazy_agg)
# Recuerda que esto es solo el cov, que se utiliza para
# reescalar luego el copy-number
# Las anotaciones en la orfs$table no han cambiado
# se cambian con esto, que está tb en el 05.filter

cons_vec <- readRDS("good_cazy_vec.rds")
orf_ids <- rownames(Arctic$orfs$table[Arctic$orfs$table$CAZy!="",])[cons_vec]
# And substitute the CAZy annotation table with blank
Arctic$orfs$table$CAZy[!(rownames(
	Arctic$orfs$table) %in% orf_ids)] <- ""

##################################
#	SAVING
##################################

saveRDS(x,file="x.rds")
save(x,y,z,file="xyz.RData")
scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/rds/file.rds /Users/alvaro/Desktop/CNB-Proyecto/Rproject/file.rds

setwd("/media/disk2/aredondo/artico022_DNA/graphs")

##################################
# SUBSETTING CONTIGS
##################################

setwd("/media/disk2/aredondo/artico022_DNA/graphs")

# Subsets
sulf.subset <- subsetFun(Arctic,"Sulfatase","PFAM")
pept.subset <- subsetORFs(Arctic,pept.orfs)
cazy.subset <- subsetORFs(Arctic,cazy.orfs)

source("analise_fun.R")

analyse.fun("sulf",sulf.subset)
analyse.fun("pept",pept.subset)
analyse.fun("cazy",cazy.subset)

##################################
# COMPARISON of USiCGs
##################################
samples <- RNAsamples
dates <- RNAdates

usiKEGG.cov <- Arctic$functions$KEGG$cov[USiKEGG,]
usiKEGG.cov <- usiKEGG.cov[,colnames(usiKEGG.cov) %in% samples]
colnames(usiKEGG.cov) <- format(as.Date(dates))

usiCOG.cov  <- Arctic$functions$COG$cov[USiCOGs,]
usiCOG.cov <- usiCOG.cov[,colnames(usiCOG.cov) %in% samples]
colnames(usiCOG.cov) <- format(as.Date(dates))

usiPFAM.cov <- Arctic$functions$PFAM$cov[USiPFAM,]
usiPFAM.cov <- usiPFAM.cov[,colnames(usiPFAM.cov) %in% samples]
colnames(usiPFAM.cov) <- format(as.Date(dates))

usicomb.cov <- interleave(t(usiKEGG.cov),t(usiCOG.cov),t(usiPFAM.cov))
# Removing outliers
usicomb.clean <- usicomb.cov
usicomb.clean[seq(3,nrow(usicomb.clean),3),c(2,6,7,13)]<-NA
# These functions are outliers in PFAM because of multiple domains

# gamma_subset$functions$COG$abund[USiCOGs,]
# Only 90, 93 and 94 have a copy in samples 9_RNA
# 90 has a more similar TPM to RecA in other samples

# Prepare for plot
usiplot <- melt(usicomb.clean)
usiplot <- cbind(usiplot,rep(c("KEGG","COG","PFAM"),nrow(usiplot)/3))
colnames(usiplot) <- c("date","gene","cov","annot")

setwd("/media/disk2/aredondo/artico022_DNA/graphs")
pdf(file="USiCGs_RNAsamples.pdf")
ggplot(usiplot, aes(x=date, y=cov, fill=annot)) +
    geom_boxplot() +
    scale_fill_manual( 
    	values = rep(c("darkgreen","aquamarine2","darkorange"),24),) +
    theme_classic() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(angle = 30,hjust=1)
    ) +
    ggtitle("USiCG coverage in RNA samples") +
    xlab("Sample") + 
    ylab("Coverage (reads per base)")
dev.off()

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/USiCGs_RNAsamples.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/USiCGs_RNAsamples.pdf

# Plot abun UScigs samples distribution







