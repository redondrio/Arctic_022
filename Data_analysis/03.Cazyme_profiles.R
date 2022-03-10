# Sections
# 1. CAZy profiles per taxon
# 2. CAZy profiles whole Arctic
# 3. Sector graphs
# 4. Beta-diversity

##################################
# CAZY PROFILES PER TAXON
##################################

colscale=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
"#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
"#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
"#8A7C64", "#599861", "#89C5FA", "#DA1724", "#74D234", "#CE51CA", "#3B4921", "#CA717C", "#CB5388", "#5B7FC7", 
"#623770", "#D5D93E", "#38A33E"
)

quart=0
btes.aggr <- aggr.cazyfam(btes_subset,"Bacteroidetes",quartile=quart)
gamma.aggr <- aggr.cazyfam(gamma_subset,"Gammaproteobacteria",quartile=quart)
#flavi.aggr <- aggr.cazyfam(flavi_subset,"Flavobacteriia",quartile=0)
alpha.aggr <- aggr.cazyfam(alpha_subset,"Alphaproteobacteria",quartile=quart)
#thiog.aggr <- aggr.cazyfam(thiog_subset,"Thioglobus",quartile=0.75)
#pelag.aggr <- aggr.cazyfam(pelag_subset,"Pelagibacterales",quartile=0.75)
#rhodo.aggr <- aggr.cazyfam(rhodo_subset,"Rhodobacterales",quartile=0.75)
archa.aggr <- aggr.cazyfam(archa_subset,"Archaea",quartile=quart)
actin.aggr <- aggr.cazyfam(actin_subset,"Actinobacteria",quartile=quart)
planc.aggr <- aggr.cazyfam(planc_subset,"Planctomycetes",quartile=quart)
verru.aggr <- aggr.cazyfam(verru_subset,"Verrucomicrobia",quartile=quart)

# Plot, including subset, function, samples and quartile threshold
pdf(file="cprof_btesGH.pdf")
ggplot(btes.aggr[btes.aggr$group=="GH",],aes(fill=tag,y=value,x=taxon)) +
	geom_bar(position="fill",stat="identity") +
	ggtitle("Taxon GH profile",subtitle=paste0("75%", " groups on DNA samples")) +
	scale_fill_manual(values=c(colscale,colscale,colscale)) +
	theme(axis.text.x = element_text(angle = 25, hjust=1))
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/cprof_btesGH.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/cprof_btesGH.pdf

for (aggregate in list(btes.aggr,gamma.aggr,flavo.aggr,thiog.aggr,
	alpha.aggr,pelag.aggr,rhodo.aggr,archa.aggr,actin.aggr,
	planc.aggr,verru.aggr)){
	n<-nrow(aggregate[aggregate$group=="GH" & aggregate$tag!="Minor",])
	print(unique(aggregate$taxon))
	print(n)
}

##################################
# CAZY PROFILES MULTIPLE TAXA
##################################

temp.aggr <- rbind(btes.aggr,gamma.aggr,alpha.aggr,archa.aggr,actin.aggr,planc.aggr,verru.aggr)

# Plot, including subset, function, samples and quartile threshold	
pdf(file=paste0("cprof_7all.pdf"))
for (fam in c("GH","GT","CE","PL")){
	print(fam)
	plot <- ggplot(temp.aggr[temp.aggr$group==fam,],aes(fill=tag,y=value,x=taxon)) +
		geom_bar(position="fill",stat="identity") +
		ggtitle(paste0("Taxon ",fam," profile"),subtitle=paste0("All families")) +
		scale_fill_manual(name="Family",values=c(colscale,colscale,colscale)) +
		labs(x="Taxon",y=paste0("% of total ",fam," copy number")) +
		theme_classic() +
		theme(axis.text.x = element_text(angle = 25, hjust=1), legend.text=element_text(size=6))
	print(plot)
}
dev.off()

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/cprof_7all.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/cprof_7all.pdf

##################################
# CAZY PROFILES WHOLE ARCTIC
##################################

# En esta parte no se hace el rowMeans porque vamos cogiendo
# una a una las muestras, por lo que el eje X va a ser el tiempo
# en lugar de los taxones como era antes

# Calculate the cazy profile per sample
arctic_aggr <- data.frame()

for (sample in DNAsamples){
	arctic_samp <- aggr.cazyfam(Arctic,"Arctic",samples=sample,quartile=0.8)
	arctic_samp$day <- sample
	arctic_aggr <- rbind(arctic_aggr,arctic_samp)
}

# Esto para que ggplot no ordene las fechas alfabéticamente
arctic_aggr$day <- factor(arctic_aggr$day, levels=unique(arctic_aggr$day))

# Plot, including subset, function, samples and quartile threshold
# Para hacer plots de todas las funciones, cambiar el fill=group
# Para hacer plots de una función, cambiar el fill=tag

pdf(file="Arctic_GH_prof.pdf")
ggplot(arctic_aggr[arctic_aggr$group=="GH",],aes(fill=tag,y=value,x=day)) +
	geom_bar(position="fill",stat="identity") + # para rescalar a 1, añadir position="fill"
	ggtitle("Arctic profile",subtitle=paste0("0.8", " groups on RNA samples")) +
	scale_fill_manual(values=c(colscale,colscale,colscale)) +
	theme(axis.text.x = element_text(angle = 25))
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/Arctic_GH_prof.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/CAZy_profiles/Arctic_GH_prof.pdf

##################################
# SECTOR GRAPHS
##################################

plots <- []
for (gp in c("GH","GT","PL","CE")){
	for (aggregate in list(btes.aggr,gamma.aggr,
		alpha.aggr,archa.aggr,actin.aggr,
		planc.aggr,verru.aggr)){
		aggregate[aggregate$tag=="Minor","tag"] <- ""
		plot <- pie(aggregate[aggregate$group=gp,"value"],
			labels = aggregate[aggregate$group==gp,"tag"])
		plots <- plots + plot
	}
}

pdf("Sector_profs.pdf")
for (plot_i in plots){print plot_i}
dev.off

# Test
pdf("Sector_profs.pdf")
btes.aggr[btes.aggr$tag=="Minor","tag"] <- ""
pie(btes.aggr[btes.aggr$group=="GH","value"],
	labels = btes.aggr[btes.aggr$group=="GH","tag"],radius=1,init.angle=90)
dev.off()

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/Sector_profs.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/Sector_profs.pdf


##################################
# CALCULATION OF DIVERSITY
##################################
library(betapart)
library(vegan)

# Shannon index (function diversity())
# Get the profiles to try and make a t-test (though it did not work)
setwd(rdsdir)
for (aggregate in list(btes.aggr,gamma.aggr,
	alpha.aggr,archa.aggr,actin.aggr,
	planc.aggr,verru.aggr)){
	write.csv(aggregate$value[aggregate$group=="GH"],file=paste0("Shannon_",aggregate$taxon[1],".csv"))
}
setwd(graphdir)
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/rds/Shannon* /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/


mmerge <- function(df,name,big){
	m <- as.matrix(df["value"])
	colnames(m) <- c(name)
	rownames(m) <- unlist(df["family"])
	merged <- merge(m,big,
		by="row.names",all.x=TRUE,all.y=TRUE)
	rownames(merged) <- merged[,1]
	return(merged[,-1])
}

for (fam in c("GH","GT","PL","CE")){
	btes.div <- as.matrix(btes.aggr[btes.aggr$group==fam,"value"])
	colnames(btes.div) <- c("Btes.")
	rownames(btes.div) <- unlist(btes.aggr[btes.aggr$group==fam,"family"])

	seven.div <- mmerge(gamma.aggr[gamma.aggr$group==fam,],"Gamma.",btes.div)
	seven.div <- mmerge(alpha.aggr[alpha.aggr$group==fam,],"Alpha.",seven.div)
	seven.div <- mmerge(actin.aggr[actin.aggr$group==fam,],"Actin.",seven.div)
	seven.div <- mmerge(archa.aggr[archa.aggr$group==fam,],"Archa.",seven.div)
	seven.div <- mmerge(planc.aggr[planc.aggr$group==fam,],"Planc.",seven.div)
	seven.div <- mmerge(verru.aggr[verru.aggr$group==fam,],"Verru.",seven.div)
	seven.div[is.na(seven.div)] <- 0

	#beta_div <- beta.pair.abund(t(seven.div), index.family = "bray")
	print(fam)
	#print(beta_div)

	mds <- metaMDS(t(seven.div))
	#print(mds)

	pdf(file=paste0("NMDS",fam,".pdf"))
	plot(mds,type="n",main=paste0("MDS respresentation of ",fam," profile distances"))
	points(mds,display="species",col="gold")
	legend(x="topright",legend=paste0(fam," families"),col="gold",pch=21)
	text(mds,display="sites",col="black")
	mtext(paste0("stress = ",round(mds$stress,2)),side=3)
	dev.off()
}

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/NMDS* /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/













