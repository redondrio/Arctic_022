##################################
# REGRESSIONS & CORRELATIONS
# Multiple functions, many taxa
##################################

# This turns off scientific notation
options(scipen = 999) # = 0 to turn back on

# Loading the metadata to see correlation with chlorophyll
metadata <- read.csv("metadata.csv")
chla <- metadata[metadata$Depth==2.5,c("DOY","Chl_a")]
# Extract the values associated with the DNA samples
chla_subset <- chla[c(1,3,3,4,6,7,8,10,13,15,16,18),2]

# Here the domains are selected, families or subfamilies
#allcazy_domains <- c(GH.domains,GT.domains,PL.domains,CE.domains)
sel_domains <- unique(fun2clean(GH.domains))

# This generates the data array with the copy numbers
cazy.trcn <- arrayTRCN(
	sub_functions=as.list(sel_domains),
	sub_levels=rep("CAZy",length(sel_domains)),
	samples=DNAsamples, #[1:12] no july [4:12] no spring
	dates=DNAdoy, #[1:12] [4:12]
	aggr_fun=sel_domains,
	gen_fun="CAZy_family",
	subsets=list(
		actin_subset,
		alpha_subset,
		gamma_subset,
		btes_subset,
		verru_subset,
		planc_subset,
		archa_subset
		),
	subset_names=c(
		"Actinobacteria",
		"Alphaproteobacteria",
		"Gammaproteobacteria",
		"Bacteroidetes",
		"Verrucomicrobia",
		"Planctomyteces",
		"Archaea"
		))

# Now the calculation of correlation
'''
# Fit a linear model to check the increase of families
# Outdated version using linear_fit
lin_fit <- function(fam){
	fam_subset <- cazy.trcn[cazy.trcn$CAZy_family==fam,]
	if (sum(fam_subset$trcn)==0){return(c(fam,rep(NA,5)))}
	linear_fit <- lm(trcn~Dates, fam_subset)
	fam_mean <- round(mean(fam_subset$trcn),5) #mean
	slope <- round(linear_fit$coefficients[2],5) #slope
	w_slope <- round(slope/fam_mean,5) #relative slope
	goodness <- round(summary(linear_fit)$r.squared,5) #R2
	pval <- round(summary(linear_fit)$coefficients[2,4],5) #pval
	return(c(fam,fam_mean,slope,w_slope,goodness,pval))
}
'''
# Function to use in mapply
# Fit a correlation model to check the relation with chlorophyll
corr_fit_2 <- function(fam,tax){
	fam_subset <- cazy.trcn[cazy.trcn$CAZy_family==fam,]
	fam_subset <- fam_subset[fam_subset$Taxon==tax,]
	fam_subset <- fam_subset[1:12,] #[1:12,] no july [4:12,] no spring
	if (sum(fam_subset$trcn)==0){
		return(c(fam,tax,rep(NA,3)))}
	corr_test <- cor.test(fam_subset[,4], log(chla_subset), method="pearson")
	fam_mean <- round(mean(fam_subset$trcn),5) #mean
	pval <- round(corr_test$p.value,5) #pval
	goodness <- round(corr_test$estimate,5) #R2
	return(c(fam,tax,fam_mean,goodness,pval))
}

# Set the subset names for the corr_fit function
subset_names=c(
		"Actinobacteria",
		"Alphaproteobacteria",
		"Gammaproteobacteria",
		"Bacteroidetes",
		"Verrucomicrobia",
		"Planctomyteces",
		"Archaea"
		)

# Mapply applies the correlationto two variables, 
# domains and taxa, and generates a data.frame
fit_list <- mapply(FUN = function(x,y) corr_fit_2(x,y), 
	x = sel_domains, y = rep(subset_names,each = length(sel_domains))) #correlation test for all families
fit_df <- as.data.frame(matrix(unlist(fit_list),nrow=length(sel_domains)*length(subset_names),byrow=T))
colnames(fit_df) <- c("CAZy","Taxon","Mean","Correlation","p_val")

# Remove blanks and order by correlation index
fit_df <- fit_df[!(is.na(levels(fit_df$Correlation)[fit_df$Correlation])),]
fit_df <- fit_df[order(as.numeric(as.character(
	fit_df$Correlation)),decreasing=TRUE),]

# Correct and threshold p-values
fit_df$cp_val <- p.adjust(
	levels(fit_df$p_val)[fit_df$p_val],
	method='fdr')
fit_df <- fit_df[fit_df$cp_val<0.05,]

# Save and load the data
write.csv2(fit_df, "CAZy_changes_all.csv", row.names = FALSE)
write.csv2(fit_df, "CAZy_changes_nojul.csv", row.names = FALSE)
write.csv2(fit_df, "CAZy_changes_centre.csv", row.names = FALSE)

fit_df <- read.csv("CAZy_changes_all.csv",sep=";") # all time points
fit_df <- read.csv("CAZy_changes_nojul.csv",sep=";") # excluding 30 Jul
fit_df <- read.csv("CAZy_changes_centre.csv",sep=";") # excluding spring as well

##################################
# HEATMAP BACTEROIDETES
##################################
cazy.trcn <- arrayTRCN(
	sub_functions=as.list(sel_domains),
	sub_levels=rep("CAZy",length(sel_domains)),
	samples=DNAsamples, #[1:12] no july [4:12] no spring
	dates=DNAdoy, #[1:12] [4:12]
	aggr_fun=sel_domains,
	gen_fun="CAZy_family",
	subsets=list(btes_subset),
	subset_names=c("Bacteroidetes"))

# Calcular correlación para ordernar
corr_fit <- function(fam){
	fam_subset <- cazy.trcn[cazy.trcn$CAZy_family==fam,]
	fam_subset <- fam_subset[1:12,] #[1:12,] no july [4:12,] no spring
	if (sum(fam_subset$trcn)==0){
		return(c(fam,rep(NA,3)))}
	corr_test <- cor.test(fam_subset[,4], chla_subset, method="spearman") #pearson before
	fam_mean <- round(mean(fam_subset$trcn),5) #mean
	pval <- round(corr_test$p.value,5) #pval
	goodness <- round(corr_test$estimate,5) #R2
	return(c(fam,fam_mean,goodness,pval))
}

fit_list <- lapply(sel_domains, FUN = function(x) corr_fit(x)) #correlation test for all families
fit_df <- as.data.frame(matrix(unlist(fit_list),nrow=length(sel_domains),byrow=T))
colnames(fit_df) <- c("CAZy","Mean","Correlation","p_val")
fit_df <- fit_df[!(is.na(levels(fit_df$Correlation)[fit_df$Correlation])),] #remove blanks
fit_df <- fit_df[order(as.numeric(as.character(
	fit_df$Correlation)),decreasing=TRUE),]

# Correct and threshold p-values
fit_df$cp_val <- p.adjust(
	levels(fit_df$p_val)[fit_df$p_val],
	method='fdr')
fit_df <- fit_df[fit_df$cp_val<0.01,]

# Sacas las CAZys ordenadas por correlación y 
# las usas otra vez para obtener un
# set vector para el heatmap 
sorted_cazy <- levels(fit_df$CAZy)[fit_df$CAZy]

cazy.trcn <- arrayTRCN(
	sub_functions=as.list(sorted_cazy),
	sub_levels=rep("CAZy",length(sorted_cazy)),
	samples=DNAsamples, #[1:12] no july [4:12] no spring
	dates=DNAdoy, #[1:12] [4:12]
	aggr_fun=sorted_cazy,
	gen_fun="CAZy_family",
	subsets=list(btes_subset),
	subset_names=c("Bacteroidetes"))

# Optional
# Normalizar sobre el máximo de cada familia
# Esto permite ver mejor las tendencias de cada familia
# pero no permite compararlas entre sí
'''for (fam in cazy.trcn$CAZy_family){
	section <- cazy.trcn[cazy.trcn$CAZy_family == fam,"trcn"]
	if (sum(section) == 0){next}
	section <- (section-min(section))/(max(section)-min(section))
    cazy.trcn[cazy.trcn$CAZy_family == fam,"trcn"] <- section
}'''

# Calcular el log para hacer la visualización más fácil
# Esto permite comparar las familias entre sí, pero no deja
# ver la tendencia de cada familia claramente
for (fam in unique(cazy.trcn$CAZy_family)){
	section <- cazy.trcn[cazy.trcn$CAZy_family == fam,"trcn"]
	if (sum(section) == 0){next}
	section <- log(section,2)
    cazy.trcn[cazy.trcn$CAZy_family == fam,"trcn"] <- section
}

# crea un factor ordenado
sorted_cazy <- factor(sorted_cazy,levels=sorted_cazy)
cazy.trcn$corr <- "neg"
for (fam in cazy.trcn$CAZy_family){
	if (fam=="GH30_3"){break} # OJO AQUÍ QUE ESTO ES MANUAL
	cazy.trcn[cazy.trcn$CAZy_family==fam,"corr"]<-"pos"
}

# Bacteroidetes plot with two colours and log copy nubmers
setwd(graphdir)
pdf(file="Btes_heatmap_Spearman.pdf")
ggplot(cazy.trcn, aes(rep(c(1,2,3,6,7,8,9,10,11,12,13,14,15),length(unique(cazy.trcn$CAZy_family))), factor(rep(sorted_cazy,each=13)), fill = corr, alpha = trcn)) + 
  geom_tile(aes(y=factor(rep(sorted_cazy,each=13)))) + 
  scale_fill_manual(values=c(neg="salmon", pos="steelblue"),name="Correlation") +
  theme_light() + theme(axis.text=element_text(size=5), legend.pos="bottom") +
  ggtitle("Evolution of GH copy-number in Bacteroidetes",subtitle="Corrected p-val (FDR) <= 0.01") +
  xlab("Dates") + ylab("Family") + labs(alpha="Log2(trcn)")
dev.off()

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/Btes_heatmap_Spearman.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/Btes_heatmap_Spearman.pdf

##################################
# HEATMAP GAMMAPROTEOBACTERIA
##################################

# Gammaproteobacteria plot with copy numbers scaled per family

cazy.trcn <- arrayTRCN(
	sub_functions=as.list(sel_domains),
	sub_levels=rep("CAZy",length(sel_domains)),
	samples=DNAsamples[1:12], #[1:12] no july [4:12] no spring
	dates=DNAdoy[1:12], #[1:12] [4:12]
	aggr_fun=sel_domains,
	gen_fun="CAZy_family",
	subsets=list(gamma_subset),
	subset_names=c("Gammaproteobacteria"))

# Optional
# Normalizar sobre el máximo de cada familia
# Esto permite ver mejor las tendencias de cada familia
# Modificado para eliminar además las filas vacías
filt_cazy <- data.frame()
for (fam in unique(cazy.trcn$CAZy_family)){
	section <- cazy.trcn[cazy.trcn$CAZy_family == fam,]
	if (sum(section[,"trcn"]) == 0){next}
	trcn <- section[,"trcn"]
	section[,"trcn"] <- (trcn-min(trcn))/(max(trcn)-min(trcn))
    filt_cazy <- rbind(filt_cazy,section)
}

setwd(graphdir)
pdf(file="Gamma_heatmap_Spearman.pdf")
# modify this to fit the    \/  number of dates \/ \/ \/
ggplot(filt_cazy, aes(rep(c(1,2,3,6,7,8,9,10,11,12,13,14),length(unique(filt_cazy$CAZy_family))), factor(filt_cazy$CAZy_family), alpha = trcn)) + 
  geom_tile(aes(y=factor(filt_cazy$CAZy_family))) + 
  theme_light() + theme(axis.text=element_text(size=5), legend.pos="bottom") +
  ggtitle("Evolution of GH copy-number in Gammaproteobacteria") +
  xlab("Dates") + ylab("Family") + labs(alpha="Max. to min. copy number")
dev.off()

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/Gamma_heatmap_Spearman.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/Gamma_heatmap_Spearman.pdf


