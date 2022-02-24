# if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
# working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/",sep="")
# setwd(working_dir)
# f<-paste(working_dir,"ABCD_release_3.0_QCed_sex_check_chrX.sexcheck",sep="")

args <- commandArgs(TRUE)
f <- args[1]
# read data
d<-read.table(f,header=TRUE)
d$SNPSEX<-as.factor(d$SNPSEX)

png(paste(f,"_hist.png",sep=""))
hist(d$F,xlab="F",main=paste("Histogram of F"))
dev.off()

library(ggplot2)
library(cowplot);theme_set(theme_cowplot())
ggplot(data=d) + geom_histogram(aes_string(x="F",fill="SNPSEX")) +
  scale_fill_brewer(palette = "Dark2",labels=c("male","female"))
ggsave(file="ABCD_release_3.0_QCed_sex_check_chrX_hist.pdf")
