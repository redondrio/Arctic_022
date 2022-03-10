# Sections
# 1. Lineplot, functions in taxon
# 2. Lineplot, function in taxa
# 3. Lineplot with chlorophyll

##################################
# TIME SERIES
# Multiple functions, one taxon
##################################

highpolar.trcn <- arrayTRCN(
	sub_functions=list(pept.domains,
		GH.domains,GT.domains,
		ref.susc),
	sub_levels=c(rep("PFAM",1),rep("CAZy",2),rep("PFAM",1)),
	samples=DNAsamples[4:13],
	dates=DNAdoy[4:13],
	aggr_fun=c("Peptidases",
		"GH","GT",
		"SusC"),
	gen_fun="Functions",
	subsets=list(polar_subset),
	subset_names=c("unclassified Polaribacter"))

lowpolar.trcn <- arrayTRCN(
	sub_functions=list("PF00884 [Sulfatase]",
		PL.domains,CE.domains,AA.domains,CBM.domains,
		ref.susd),
	sub_levels=c(rep("PFAM",1),rep("CAZy",4),rep("PFAM",1)),
	samples=DNAsamples[4:13],
	dates=DNAdoy[4:13],
	aggr_fun=c("Sulfatases",
		"PL","CE","AA","CBM",
		"SusD"),
	gen_fun="Functions",
	subsets=list(polar_subset),
	subset_names=c("unclassified Polaribacter"))

setwd("/media/disk2/aredondo/artico022_DNA/graphs")
pdf(file="time_trcn_polar.pdf")
high_plot <- ggplot(highpolar.trcn, aes(x = Dates, y = trcn, color=Functions)) + 
				geom_smooth(method=lm) + geom_point() + theme_minimal() +
				labs(title = "Function abundance in unclassified Polaribacter",
					y = "Taxon-rescaled copy number",
					x = "Day of year")
print(high_plot)
low_plot <- ggplot(lowpolar.trcn, aes(x = Dates, y = trcn, color=Functions)) + 
				geom_smooth(method=lm) + geom_point() + theme_minimal() +
				labs(title = "Function abundance in unclassified Polaribacter",
					y = "Taxon-rescaled copy number",
					x = "Day of year")
print(low_plot)
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/time_trcn_polar.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/time_trcn_polar.pdf


##################################
# TIME SERIES
# One function, multiple taxons
##################################

allcazy_domains <- c(GH.domains,GT.domains,PL.domains,CE.domains)
sel_domains <- unique(fun2fam(GH.domains))
cazy.trcn <- arrayTRCN(
	sub_functions=as.list(sel_domains),
	sub_levels=rep("CAZy",length(sel_domains)),
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

setwd("/media/disk2/aredondo/artico022_DNA/graphs")
pdf(file="Trcn7_per_GH.pdf")
for (cazy_fam in sort(sel_domains)){
	if (sum(cazy.trcn[cazy.trcn$CAZy_family==cazy_fam,"trcn"])==0){next}
	cazy_plot <- ggplot(cazy.trcn[cazy.trcn$CAZy_family==cazy_fam,],
	 aes(x = Dates, y = trcn, color=Taxon)) + 
					#geom_smooth(method=lm,se=FALSE) + 
					geom_line() + geom_point() + theme_minimal() +
					labs(title = paste0(cazy_fam," trcn in time"),
						y = "Taxon-rescaled copy number",
						x = "Day of year")
	print(cazy_plot)
	print(cazy_fam)
}
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/Trcn7_per_GH.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/Trcn7_per_GH.pdf

##################################
# TIME SERIES
# Taxa with chlorophyll background
##################################

# This turns off scientific notation
options(scipen = 999) # = 0 to turn back on

# Loading the metadata to see correlation with chlorophyll
metadata <- read.csv("metadata.csv")
chla <- metadata[metadata$Depth==2.5,c("DOY","Chl_a")]
# Extract the values associated with the dates of DNA samples
chla_subset <- as.data.frame(
	chla[c(1,3,4,6,7,8,10,13,15,16,18),])

# Extract the percentages of taxa
sk_plot <- plotTaxonomy(Arctic,
	rank = "superkingdom",
	count = "percent",
	tax = "Archaea",
	other = F,
	samples = DNAsamples)

ph_plot <- plotTaxonomy(Arctic,
	rank = "phylum",
	count = "percent",
	tax = c("Bacteroidetes","Planctomycetes",
		"Verrucomicrobia","Actinobacteria"),
	other = F,
	samples = DNAsamples)

cl_plot <- plotTaxonomy(Arctic,
	rank = "class",
	count = "percent",
	tax = c("Alphaproteobacteria",
		"Gammaproteobacteria"),
	other = F,
	samples = DNAsamples)

or_plot <- plotTaxonomy(Arctic,
	rank = "order",
	count = "percent",
	tax = c("Flavobacteriales","Pelagibacterales",
		"Rhodobacterales","Cellvibrionales"),
	other = F,
	samples = DNAsamples)

gn_plot <- plotTaxonomy(Arctic,
	rank = "genus",
	count = "percent",
	tax = c("SAR86.cluster..no.genus.in.NCBI.",
		"Candidatus.Pelagibacter",
		"Candidatus.Thioglobus","Polaribacter"),
	other = F,
	samples = DNAsamples)

# Colours
'''
New consistent colours
Btes      Archa     Alpha     Gamma     Planct    Verru     Actin
"#CD672F","#E01853","#1F74CD","#5D478B","#F6BFCA","#8B461D","#FCD84A"
'''

ph_colour <- c("#CD672F", #Btes
			 "#F6BFCA", #Planct
			 "#8B461D", #Verru
			 "#FCD84A", #Actin
			 "#E01853", #Archa
			 "#1F74CD", #Alpha
			 "#5D478B"  #Gamma
			 )

# Prepare data for plotting
plot1_data <- rbind(ph_plot$data,sk_plot$data)
plot1_data <- rbind(plot1_data,cl_plot$data)
plot1_data$DOY <- DNAdoy
plot1_data <- subset(plot1_data,sample!="30_Jul")

plot2_data <- rbind(or_plot$data)
plot2_data$DOY <- DNAdoy
plot2_data <- subset(plot2_data,sample!="30_Jul")

plot3_data <- rbind(gn_plot$data)
plot3_data$DOY <- DNAdoy
plot3_data <- subset(plot3_data,sample!="30_Jul")
plot3_data$item <- c(rep("SAR86 cluster",12),
					 rep("C. Pelagibacter",12),
					 rep("C. Thioglobus",12),
					 rep("Polaribacter",12))

# Plot
plot1 <- ggplot(plot1_data,aes(x=DOY)) + geom_area(
	data=chla_subset, aes(
		x=DOY, y=Chl_a*max(plot1_data$abun)/max(chla_subset$Chl_a)),
		fill="#b3dec1", colour="gray") + scale_y_continuous(
	sec.axis = sec_axis(
		~./max(plot1_data$abun)*max(chla_subset$Chl_a), name="Chlorophyll A")) + geom_line(
	data=plot1_data, aes(
		x=DOY, y=abun, group=item, colour=item)) + geom_point(
	data=plot1_data, aes(
		x=DOY, y=abun, colour=item)
		) + theme_classic() + scale_color_brewer(palette="Dark2") + labs(
		x="Day of year",y="Abundance (%of reads)", colour="Clade") + theme(legend.position="top")

plot2 <- ggplot(plot2_data,aes(x=DOY)) + geom_area(
	data=chla_subset, aes(
		x=DOY, y=Chl_a*max(plot2_data$abun)/max(chla_subset$Chl_a)),
	fill="#b3dec1", colour="gray") + scale_y_continuous(
	sec.axis = sec_axis(
		~./max(plot2_data$abun)*max(chla_subset$Chl_a), name="Chlorophyll A")) + geom_line(
	data=plot2_data, aes(
		x=DOY, y=abun, group=item, colour=item)) + geom_point(
	data=plot2_data, aes(
		x=DOY, y=abun, colour=item)
		) + theme_classic() + scale_color_brewer(palette="Dark2") + labs(
		x="Day of year",y="Abundance (%of reads)", colour="Clade") + theme(legend.position="top")

plot3 <- ggplot(plot3_data,aes(x=DOY)) + geom_area(
	data=chla_subset, aes(
		x=DOY, y=Chl_a*max(plot3_data$abun)/max(chla_subset$Chl_a)),
	fill="#b3dec1", colour="gray") + scale_y_continuous(
	sec.axis = sec_axis(
		~./max(plot3_data$abun)*max(chla_subset$Chl_a), name="Chlorophyll A")) + geom_line(
	data=plot3_data, aes(
		x=DOY, y=abun, group=item, colour=item)) + geom_point(
	data=plot3_data, aes(
		x=DOY, y=abun, colour=item)
		) + theme_classic() + scale_color_brewer(palette="Dark2") + labs(
		x="Day of year",y="Abundance (%of reads)", colour="Clade") + theme(legend.position="top")

pdf(file="Teelings_small1.pdf",width=5,height=8)
ggplot(plot1_data,aes(x=DOY)) + 
	geom_area(
	data=chla_subset, aes(
		x=DOY, y=Chl_a*max(plot1_data$abun)/max(chla_subset$Chl_a)),
		fill="#C0F1D3", colour="gray") + 
	scale_y_continuous(
	sec.axis = sec_axis(
		~./max(plot1_data$abun)*max(chla_subset$Chl_a), name="Chlorophyll A")) + 
	geom_line(
	data=plot1_data, aes(
		x=DOY, y=abun, group=item, colour=item), size=0.6) + 
	geom_point(
	data=plot1_data, aes(
		x=DOY, y=abun, colour=item), size=2) + 
	theme_classic() + scale_color_manual(values=ph_colour) + 
	labs(
		x="Day of year",y="Abundance (% of reads)", colour="Clade") + 
	theme(legend.position="top") + 
	guides(colour = guide_legend(ncol = 2))

#print(plot1)
#print(plot2)
#print(plot3)
dev.off()

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/Teelings_small1.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/Teelings_small1.pdf


















