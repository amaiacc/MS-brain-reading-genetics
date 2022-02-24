# library(ggpubr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(tidyr)
#----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
# define pattern for this run
root="ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned_rm05_adj"

# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/GCTA/",sep="")
reml_dir=paste(working_dir,"/reml/",sep="")
out_dir=paste(working_dir,"/output/",sep="")
#
dir.create(file.path(out_dir), showWarnings = FALSE)
setwd(out_dir)
#----------------------------------------------------------------------
h2<-read.table(paste(reml_dir,"hsq_h2_summary_",root,".table",sep=""),stringsAsFactors = FALSE)
colnames(h2)<-c("file","n","h2","se","pval")
h2$file<-gsub(":n","",h2$file)
h2$p1<-gsub(".hsq","",h2$file) %>% gsub("_residuals|_covs","",.) 
h2$stat<-"h2"
h2$sample<-strsplit(root,"QCed_") %>% sapply("[[",2)
h2$run<-NA
h2$run[grep("residuals",h2$file)]<-"residuals"
h2$run[grep("covs",h2$file)]<-"covs"
h2$adj<-"-"
h2$adj[grep("adjGlobal",h2$file)]<-"adjGlobal"

# define categories for the phenotypes
h2$category<-NA
h2$category[grep("smri*.*total",h2$p1)]<-"Brain global"
h2$category[grep("smri",h2$p1)[grep("total",h2$p1[grep("smri",h2$p1)],invert=TRUE)]]<-"Brain ROI"
h2$category[grep("smri",h2$p1,invert=TRUE)]<-"Cognitive"
# for brain regions: measure type, hemisphere and region name
h2$measure<-h2$region<-h2$hemisphere<-NA
w<-which(h2$category!="Cognitive")
h2$hemisphere[w][grep("lh|LH",h2$p1[w])]<-"L"
h2$hemisphere[w][grep("rh|RH",h2$p1[w])]<-"R"
h2$region[w]<-gsub(".lh|.rh","",h2$p1[w])
h2$measure[w][grep("area",h2$p1[w])]<-"AREA"
h2$measure[w][grep("thick",h2$p1[w])]<-"THICKNESS"
# save 
write.csv(h2,  file=paste(out_dir,"gcta_estimates_h2_",root,".csv",sep=""),row.names = FALSE,quote=TRUE)
#-------------------------
rg<-read.table(paste(reml_dir,"hsq_rg_summary_",root,".table",sep=""),stringsAsFactors = FALSE)
colnames(rg)<-c("file","rg","se","pval")
rg$file<-gsub(":rG","",rg$file)
rg$p1<-gsub(".hsq","",rg$file) %>% strsplit(.,"_residuals|_covs") %>% sapply("[[",1) %>% gsub("^_|_$","",.) %>% gsub("adjGlobal_","",.)
rg$p2<-gsub(".hsq","",rg$file) %>% strsplit(.,"_residuals|_covs") %>% sapply("[[",2) %>% gsub("^_|_$","",.) %>% gsub("adjGlobal_","",.)
rg$stat<-"rg"
rg$sample<-strsplit(root,"QCed_") %>% sapply("[[",2)
rg$run<-NA
rg$run[grep("residuals",rg$file)]<-"residuals"
rg$run[grep("covs",rg$file)]<-"covs"
rg$adj<-"-"
rg$adj[grep("adjGlobal",rg$file)]<-"adjGlobal"

# order pheno cols, so that we can remove duplicates
rg[,c("p1","p2")]<-t(apply(rg[,c("p1","p2")], 1, sort))
# delete duplicated
w<-which(duplicated(rg[,c("rg","se","p1","p2")]))

rg[w,]
if (length(w)>0){rg<-rg[-w,]}
rm(w)

#
write.csv(rg,  file=paste(out_dir,"gcta_estimates_rg_",root,".csv",sep=""),row.names = FALSE,quote=TRUE)

