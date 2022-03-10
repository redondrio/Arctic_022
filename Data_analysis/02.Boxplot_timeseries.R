##################################
# BOXPLOTS
##################################

# Array with the trcn of the functions selected in the taxons selected
# Takes the functions you want in the X axis as lists of domains
# that are later aggregated with a group name
# Takes the subsets you want as different boxes
# and the names for the legend
# Returns a dataframe that can be directly plotted with ggplot	

#Prueba a hacerlas con todas las muestas, no solo ADN
fun.trcn <- arrayTRCN(
	sub_functions=list(pept.domains,"PF00884 [Sulfatase]",AA.domains,CBM.domains),
	sub_levels=c("PFAM","PFAM","CAZy","CAZy"),
	aggr_fun=c("Peptidases","Sulfatases","AA","CBM"),
	gen_fun="Functions",
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

cazy.trcn <- arrayTRCN(
	sub_functions=list(GH.domains,GT.domains,PL.domains,CE.domains),
	sub_levels=rep("CAZy",4),
	aggr_fun=c("GH","GT","PL","CE"),
	gen_fun="CAZy_group",
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

sus.trcn <- arrayTRCN(
	sub_functions=list(ref.susd,ref.susc),
	sub_levels=rep("PFAM",2),
	aggr_fun=c("SusD homologues","SusC homologues"),
	gen_fun="PUL_elements",
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
sus.trcn$Taxon <- factor(sus.trcn$Taxon,levels=c("Actinobacteria",
													"Archaea",
													"Alphaproteobacteria",
													"Gammaproteobacteria",
													"Bacteroidetes",
													"Planctomyteces",
													"Verrucomicrobia")) 


wide.trcn <- arrayTRCN(
	sub_functions=list(GH.domains,GT.domains,PL.domains,CE.domains,pept.domains,"PF00884 [Sulfatase]"),
	sub_levels=c(rep("CAZy",4),"PFAM","PFAM"),
	aggr_fun=c("GH","GT","PL","CE","Peptidases","Sulfatases"),
	gen_fun="Functions",
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

wide.trcn$Functions <- factor(wide.trcn$Functions,levels=c("GH","GT","PL","CE","Peptidases","Sulfatases")) 
wide.trcn$Taxon <- factor(wide.trcn$Taxon,levels=c("Actinobacteria",
													"Archaea",
													"Alphaproteobacteria",
													"Gammaproteobacteria",
													"Bacteroidetes",
													"Planctomyteces",
													"Verrucomicrobia")) 

#aggregate(fun.trcn,by=fun.trcn[,c(2,3)],FUN=mean)

#library(reshape)
#library(ggplot2)

##################################
# BOXPLOTS
##################################

'''
New consistent colours
Btes      Archa     Alpha     Gamma     Planct    Verru     Actin
"#CD672F","#E01853","#1F74CD","#5D478B","#F6BFCA","#8B461D","#FCD84A"
'''

az_colour <- c("#FCD84A", #Actin
			"#E01853", #Archa
			"#1F74CD", #Alpha
			"#5D478B", #Gamma
			"#CD672F", #Btes
			"#F6BFCA", #Planc
			"#8B461D"  #Verru
			 )

setwd("/media/disk2/aredondo/artico022_DNA/graphs")
pdf(file="boxplots7_fun.pdf")
ggplot(fun.trcn, aes(x = Functions, y = trcn, color=Taxon)) + 
geom_boxplot() + theme_light() + theme(legend.position="bottom") +
ylab("Copy-number")
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/boxplots7_fun.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/boxplots7_fun.pdf

setwd("/media/disk2/aredondo/artico022_DNA/graphs")
pdf(file="boxplots7_cazy.pdf")
ggplot(cazy.trcn, aes(x = CAZy_group, y = trcn, color=Taxon)) + 
geom_boxplot() + ylim(0,50) + theme_light() + theme(legend.position="bottom") + 
xlab("CAZy group") + ylab("Copy-number")
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/boxplots7_cazy.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/boxplots7_cazy.pdf

setwd("/media/disk2/aredondo/artico022_DNA/graphs")
pdf(file="boxplots7_sus.pdf",width=4,height=6)
ggplot(sus.trcn, aes(x = PUL_elements, y = trcn, color=Taxon, fill=Taxon)) + 
geom_boxplot(alpha=0.2) + theme_light() + theme(legend.position="bottom") + 
scale_color_manual(values = az_colour) +
scale_fill_manual(values = az_colour) +
xlab("PUL element") + ylab("Copy-number") +
guides(colour = guide_legend(ncol = 2))
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/boxplots7_sus.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/boxplots7_sus.pdf

setwd("/media/disk2/aredondo/artico022_DNA/graphs")
pdf(file="boxplots7_wide.pdf",width=8,height=4)
ggplot(wide.trcn, aes(x = Functions, y = trcn, color=Taxon, fill=Taxon)) + 
geom_boxplot(alpha=0.2) + theme_light() + theme(legend.position="bottom") + 
scale_color_manual(values = az_colour) +
scale_fill_manual(values = az_colour) +
xlab("Function") + ylab("Copy-number")
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/boxplots7_wide.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/boxplots7_wide.pdf

