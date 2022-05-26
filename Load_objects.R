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
library(vegan)
tf <- readRDS("tandem_frame.rds")
susc.allindx <- readRDS("susc_allindx.rds")
susd.allindx <- readRDS("susd_allindx.rds")
rdsdir <- "/media/disk2/aredondo/artico022_DNA/rds"
graphdir <- "/media/disk2/aredondo/artico022_DNA/graphs"

# Esto es para modificar el objeto Arctic con
# el resultado de filtrar las cazymes
# Recuerda que esto es solo el cov, que se utiliza para
# reescalar luego el copy-number
cazy_agg <- readRDS("filt_cazy_cov.rds")
Arctic$functions$CAZy$cov <- as.matrix(cazy_agg)
# Las anotaciones en la orfs$table no han cambiado
# se cambian con esto, que está tb en el 05.filter
cons_vec <- readRDS("good_cazy_vec.rds")
orf_ids <- rownames(Arctic$orfs$table[Arctic$orfs$table$CAZy!="",])[cons_vec]
# And substitute the CAZy annotation table with blank
Arctic$orfs$table$CAZy[!(rownames(
	Arctic$orfs$table) %in% orf_ids)] <- ""
