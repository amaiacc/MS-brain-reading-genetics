library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
options(stringsAsFactors = FALSE)
# 

args <- commandArgs(TRUE)
root <- args[1]
pcx <- as.numeric(args[2])
pcy <- as.numeric(args[3])

PCAEVEC<-read.table(paste(root,".pca.evec",sep=""),head=F)
# get number of pcs from file
n=NCOL(PCAEVEC)-3
#
colnames(PCAEVEC)<-c("FID","IID",
                     paste("PC",1:n,sep=""),
                     "Pheno")

# read info about batches if available, if not avail, it will be the same as r
r<-gsub(".pop_strat|.pop_strat_outliers","",root)
if (file.exists(pattern=paste(r,".batches",sep=""))){
  batch<-read.table(paste(r,".batches",sep=""),header=TRUE)
  PCAEVEC<-merge(PCAEVEC,batch,by=c("FID","IID"),all.x=TRUE)
  PCAEVEC$Pheno<-PCAEVEC$batch
} else {
  PCAEVEC$Pheno<-r
  }
rm(r)

# plot
pcx<-2+pcx
pcy<-2+pcy
# pdf(paste(root,"_PC",(pcx-2),"_PC",(pcy-2),".pdf",sep=""))
png(paste(root,"_PC",(pcx-2),"_PC",(pcy-2),".png",sep=""))
print(
  qplot(PCAEVEC[,pcx],PCAEVEC[,pcy], data=PCAEVEC, color=Pheno) + xlab(paste("PC",(pcx-2),sep="")) + ylab(paste("PC",(pcy-2),sep="")) + theme(legend.position="bottom")
)
dev.off()

# ## loop
# for (x in 1:3){
#     pcx=x
#     pcy=x+1
#     print(paste(pcx, "vs",pcy))
#     pcx<-2+pcx
#     pcy<-2+pcy
#     pdf(paste(root,"_PC",(pcx-2),"_PC",(pcy-2),".pdf",sep=""))
#     print(
#      qplot(PCAEVEC[,pcx],PCAEVEC[,pcy], data=PCAEVEC, color=Pheno) + xlab(paste("PC",(pcx-2),sep="")) + ylab(paste("PC",(pcy-2),sep=""))
#     )
#     dev.off()
# }

