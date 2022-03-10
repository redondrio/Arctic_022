####################################
# CHECK THE EVOLUTION OF CAZYME
# ANNOTATIONS WITH VARYING
# THRESHOLDS AND PLOT
####################################

# R now
library(ggplot2)
library(reshape2)
cazyfam <- "GH"
genomeid_list <- c("153266.4",
"269798.16",
"402612.5",
"411154.11",
"1347342.6",
"376686.10",
"1336803.3",
"1250005.3",
"1836467.10",
"248905.6",
"248906.8",
"313598.6",
"1977054.20",
"641691.6",
"1897614.9",
"2268138.3",
"1897614.10",
"313590.8")

eval_df <- data.frame()
id_df <- data.frame()
cov_df <- data.frame()

for (genomeid in genomeid_list){
	setwd(paste0("/media/disk2/aredondo/artico022_DNA/genomes/",genomeid,"/results/GH_quality"))
	# Load the diamond results and the reference annotations
	# for GTs as extracted before
	dmnd.res <- read.table(file="dmnd_GH.tsv")
	ref.ids <- read.table(file="ref_GH_ids.txt")
	names(dmnd.res) <- c("query","query_length","hit","hit_length",
		"identity","matches","e-value","score",
		"query_start","query_end","hit_start","hit_end")
	dmnd.res$cov <- (dmnd.res$query_end-dmnd.res$query_start)/dmnd.res$query_length
	# Let's see how different thresholds change the number
	# of recovered GTs
	# This generates two vectors, rec and acc, for the recall
	# and precision of annotations with different id thresholds
	
	# This one is for E-value
	thresholds_eval <- as.numeric(seq(10,190,10))
	rec_eval <- c()
	acc_eval <- c()
	for (h in thresholds_eval){
		sqm_pos <- (length(unique(dmnd.res[dmnd.res[,7]<10**(-h),1])))
		true_pos <- (sum(ref.ids[,1] %in% dmnd.res[dmnd.res[,7]<10**(-h),1]))
		recovery <- round(true_pos/nrow(ref.ids),2)
		accuracy <- round(true_pos/sqm_pos,2)
		rec_eval <- c(rec_eval,recovery)
		acc_eval <- c(acc_eval,accuracy)
	}

	# This one is for ID%
	thresholds_id <- as.numeric(seq(30,100,2))
	rec_id <- c()
	acc_id <- c()
	for (h in thresholds_id){
		sqm_pos <- (length(unique(dmnd.res[dmnd.res[,5]>=h,1])))
		true_pos <- (sum(ref.ids[,1] %in% dmnd.res[dmnd.res[,5]>=h,1]))
		recovery <- round(true_pos/nrow(ref.ids),2)
		accuracy <- round(true_pos/sqm_pos,2)
		rec_id <- c(rec_id,recovery)
		acc_id <- c(acc_id,accuracy)
	}

	# This one is for coverage
	thresholds_cov <- as.numeric(seq(0,1.01,0.02))
	rec_cov <- c()
	acc_cov <- c()
	for (h in thresholds_cov){
		sqm_pos <- (length(unique(dmnd.res[dmnd.res[,13]>=h,1])))
		true_pos <- (sum(ref.ids[,1] %in% dmnd.res[dmnd.res[,13]>=h,1]))
		recovery <- round(true_pos/nrow(ref.ids),2)
		accuracy <- round(true_pos/sqm_pos,2)
		rec_cov <- c(rec_cov,recovery)
		acc_cov <- c(acc_cov,accuracy)
	}


	# This will turn them into dfs
	data_eval <- data.frame(
						id = genomeid,
						x = thresholds_eval,
	                    Recall = rec_eval,
	                    Precision = acc_eval)
	# This will turn them into dfs
	data_id <- data.frame(
						id = genomeid,
						x = thresholds_id,
	                    Recall = rec_id,
	                    Precision = acc_id)
	# This will turn them into dfs
	data_cov <- data.frame(
						id = genomeid,
						x = thresholds_cov,
	                    Recall = rec_cov,
	                    Precision = acc_cov)

	# Add to the general df
	eval_df <- rbind(eval_df,data_eval)
	id_df <- rbind(id_df,data_id)
	cov_df <- rbind(cov_df,data_cov)

}

# Prepare them for ggplot
melted_eval <- melt(eval_df,measure.vars=c("Recall","Precision"))
melted_eval$key <- paste0(melted_eval$id,melted_eval$variable)

melted_id <- melt(id_df,measure.vars=c("Recall","Precision"))
melted_id$key <- paste0(melted_id$id,melted_id$variable)

melted_cov <- melt(cov_df,measure.vars=c("Recall","Precision"))
melted_cov$key <- paste0(melted_cov$id,melted_cov$variable)

setwd(paste0("/media/disk2/aredondo/artico022_DNA/genomes"))
pdf(file="GH_quality.pdf")

ggplot(melted_eval[melted_eval$variable=="Recall",],aes(x=x,y=value)) +
	geom_line(aes(linetype = variable,color=id)) +
	geom_smooth(method=NULL, color="blue", se=TRUE) +
	labs(y="% recovered", x = "E-value threshold") +
	ylim(0,1)
ggplot(melted_eval[melted_eval$variable=="Precision",],aes(x=x,y=value)) +
	geom_line(aes(linetype = variable,color=id)) +
	geom_smooth(method=NULL, color="blue", se=TRUE) +
	labs(y="% precision", x = "E-value threshold") +
	ylim(0,1)

ggplot(melted_id[melted_id$variable=="Recall",],aes(x=x,y=value)) +
	geom_line(aes(linetype = variable,color=id)) +
	geom_smooth(method=NULL, color="blue", se=TRUE) +
	labs(y="% recovered", x = "ID% threshold") +
	ylim(0,1)
ggplot(melted_id[melted_id$variable=="Precision",],aes(x=x,y=value)) +
	geom_line(aes(linetype = variable,color=id)) +
	geom_smooth(method=NULL, color="blue", se=TRUE) +
	labs(y="% precision", x = "ID% threshold") +
	ylim(0,1)

ggplot(melted_cov[melted_cov$variable=="Recall",],aes(x=x,y=value)) +
	geom_line(aes(linetype = variable,color=id)) +
	geom_smooth(method=NULL, color="blue", se=TRUE) +
	labs(y="% recovered", x = "Coverage threshold") +
	ylim(0,1)
ggplot(melted_cov[melted_cov$variable=="Precision",],aes(x=x,y=value)) +
	geom_line(aes(linetype = variable,color=id)) +
	geom_smooth(method=NULL, color="blue", se=TRUE) +
	labs(y="% precision", x = "Coverage threshold") +
	ylim(0,1)
dev.off()
#scp aredondo@zobel.cnb.csic.es:/media/disk2/aredondo/artico022_DNA/genomes/GH_quality.pdf /Users/alvaro/Desktop/CNB-Proyecto/Rproject/Graphs/GH_quality.pdf



