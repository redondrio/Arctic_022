##################################
# RESCALING TO TRCN
##################################

##### big_subset<-subsetTax(Arctic,"phylum","Bacteroidetes",rescale_tpm=F,rescale_copy_number=F)

# Careful here, rescaling DNA samples gives you the number of gene
# copies per genome, but rescaling RNA samples gives you the expression
# ratio related to single copy genes, which doesn't mean much

setwd(rdsdir)

# Load saved objects
btes_subset <- readRDS(file="btes_sublite_norm.rds")
gamma_subset <- readRDS(file="gamma_sublite_norm.rds")
alpha_subset <- readRDS(file="alpha_sublite_norm.rds")
archa_subset <- readRDS(file="archa_sublite_norm.rds")
actin_subset <- readRDS(file="actin_sublite_norm.rds")
planc_subset <- readRDS(file="planc_sublite_norm.rds")
verru_subset <- readRDS(file="verru_sublite_norm.rds")

# Already created the subsets, just leave it here as a reminder

'''
# Create, rescale and save subsets of a given taxon
# Re-did the whole thing with rescaled TPMs
btes_subset  <- save.litesub(Arctic,"btes","phylum","Bacteroidetes")
flavo_subset <- save.litesub(Arctic,"flavo","order","Flavobacteriales")
thiog_subset <- save.litesub(Arctic,"thiog","species","Candidatus Thioglobus")
gamma_subset <- save.litesub(Arctic,"gamma","class","Gammaproteobacteria")
pelag_subset <- save.litesub(Arctic,"pelag","order","Pelagibacterales")
rhodo_subset <- save.litesub(Arctic,"rhodo","order","Rhodobacterales")
alpha_subset <- save.litesub(Arctic,"alpha","class","Alphaproteobacteria")
actin_subset <- save.litesub(Arctic,"actin","class","Actinobacteria")
archa_subset <- save.litesub(Arctic,"archa","superkingdom","Archaea")
planc_subset <- save.litesub(Arctic,"planc","phylum","Planctomycetes")
verru_subset <- save.litesub(Arctic,"verru","phylum","Verrucomicrobia")
flavi_subset <- save.litesub(Arctic,"flavi","class","Flavobacteriia")
flavu_subset <- save.litesub(Arctic,"flavu","genus","Unclassified Flavobacteriaceae")
flave_subset <- save.litesub(Arctic,"flave","family","Flavobacteriaceae")
polar_subset <- save.litesub(Arctic,"polar","genus","Polaribacter")
#useless because there are too few reads
#cellv_subset <- save.litesub(gamma_subset,"cellv","order","Cellvibrionales")
'''