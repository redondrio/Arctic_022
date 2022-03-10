##################################
# PIECHARTS
##################################

# Use an array with the trcn of the functions selected in the taxons selected
# Generates a grid of piecharts showing the copy number of each taxon
# in the functions you want in the X axis as lists of CAZy domains
# and the samples you want in the Y axis

pieplot_domains <- c("GH9","GH16","GH17","GH30","GH55","GH81","GH128") #laminarin
pieplot_domains <- c("GH26","GH38","GH76","GH92","GH99","GH113","GH125", "GH148") #mannanase
pieplot_domains <- c("GH5","GH6","GH7","GH8","GH9","GH12","GH44","GH45","GH48","GH124") #cellulase
pieplot_domains <- c("GH5","GH13","GH15","GH57","GH97","GH119","GH126") #starch
pieplot_domains <- c("GH27","GH35","GH36","GH42", "GH29","GH95") #galactanase, fucosidase
pieplot_domains <- c("GH2","GH29","GH95","GH42","GH3","GH16","GH92","GH17","GH30","GH43",
					 "GH1","GH31","GH5","GH13","GH20","GH81","GH27","GH65") #Teelings

# Names to write on the first column
pieplot_names <- c(
	GH1 = "b-Glycosidase",
	GH2 = "Gal/GlcA",
	GH3 = "Glycosidase",
	GH5 = "Glc/Man",
	GH6 = "Cellulase",
	GH13 = "a-Amylase",
	GH16 = "Laminarinase",
	GH17 = "Laminarinase",
	GH20 = "Chitinase",
	GH27 = "a-Galactosidase",
	GH29 = "Fucosidase",
	GH30 = "b-Glucosidase",
	GH31 = "a-Glycosidase",
	GH42 = "b-Galactosidase",
	GH43 = "a-L-Arabinases",
	GH65 = "Phosphorylase",
	GH81 = "Laminarinase",
	GH92 = "a-Mannosidase",
	GH95 = "Fucosidase"
	)

# Samples as columns
samples <- DNAdoy

# Taxa in piecharts
taxa <- c(
		"Actinobacteria",
		"Alphaproteobacteria",
		"Gammaproteobacteria",
		"Bacteroidetes",
		"Verrucomicrobia",
		"Planctomyteces",
		"Archaea"
		)

pie.trcn <- arrayTRCN(
	sub_functions=pieplot_domains,
	sub_levels=rep("CAZy",length(pieplot_domains)),
	aggr_fun=pieplot_domains,
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
	subset_names=taxa)

# GRID OF PIECHARTS
l_dom <- length(pieplot_domains) # number of rows
l_sam <- length(samples) # number of columns

# Colours as a named vector
collist <- c("#FCD84A", #Actin
			 "#1F74CD", #Alpha
			 "#5D478B", #Gamma
			 "#CD672F", #Btes
			 "#8B461D", #Verru
			 "#F6BFCA", #Planc
			 "#E01853"  #Archa
			 )
names(collist) <- taxa

'''
New consistent colours
Btes      Archa     Alpha     Gamma     Planct    Verru     Actin
"#CD672F","#E01853","#1F74CD","#5D478B","#F6BFCA","#8B461D","#FCD84A"
'''

# PLOT

pdf("piegrid_cazyme_scaled.pdf")

layout_matrix <- matrix(1:((l_dom+1)*(l_sam+1)),ncol=(l_sam+1),byrow=TRUE) # add one row for the names
layout(layout_matrix) # activate layout

par(mar=c(0,0,0,0))
for (cazyme in pieplot_domains){ # for each row
	plot(0,bty="n",type="n",xaxt="n",yaxt="n",ann=F) # empty plot
	text(1,0,paste0(
		as.character(cazyme),"\n",pieplot_names[cazyme]
		),col="black",cex=0.5) # print rowname
	max_cn <- max(Arctic$functions$CAZy$copy_number[cazyme,DNAsamples]) # to scale radius per cazyme
	min_cn <- min(Arctic$functions$CAZy$copy_number[cazyme,DNAsamples])
  	for (sample_i in seq(1:length(samples))){ # for each column
  		data <- pie.trcn[which(pie.trcn$Functions==cazyme & pie.trcn$Dates==samples[sample_i]),4]
  		#max_cn <- max(Arctic$functions$CAZy$copy_number[pieplot_domains,DNAsamples[sample_i]]) # to scale radius per sample
  		#min_cn <- min(Arctic$functions$CAZy$copy_number[pieplot_domains,DNAsamples[sample_i]])
  		ex_taxa <- taxa[data!=0] # select taxa with >0 trcn
  		if (length(ex_taxa)==0){
  			plot(0,bty="n",type="n",xaxt="n",yaxt="n",ann=F) # empty plot if no data
  		}
  		else{
	  		data <- data[data!=0] # select data >0 to avoid small sectors in the pie
	  		rad <- (Arctic$functions$CAZy$copy_number[cazyme,DNAsamples[sample_i]]-min_cn)/(max_cn-min_cn)*(1-0.2)+0.2 # to have a minimum radius
			pie(data,label=NA,col=collist[ex_taxa],lty="blank",radius=rad)
			# pie the copy numbers of the selected cazyme and sample
		}
	}
}

plot(0,bty="n",type="n",xaxt="n",yaxt="n",ann=F) # empty plot for the lower left corner

for (sample in DNAsamples){ # print column names
	plot(0,bty="n",type="n",xaxt="n",yaxt="n",ann=F) # empty plot
	text(1,0,as.character(sample),col="black") # print colname
}

dev.off()

#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/graphs/piegrid_cazyme_scaled.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/piegrid_cazyme_scaled.pdf











