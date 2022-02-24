#!/opt/R-4.0.3/bin/R

# clean workspace
rm(list=ls())

# parameters
# args <- commandArgs(TRUE)

# define libraries
library(dplyr); library(tidyr)

# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
pgeno_dir=paste(dir,"resources/datasets/ABCD/data/primary_data/ABCDgenomicsDataV30/genomics_sample03/ABCD_genotype/",sep="")
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/",sep="")
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",sep="")
base_dir=paste(geno_dir,"/GCTA/pheno/",sep="")
if (!dir.exists(base_dir)){dir.create(base_dir)}

setwd(pheno_dir)

# read data
info<-read.table(paste(pgeno_dir,"ABCD_release3.0_.batch_info.txt",sep=""),header=TRUE)
dat <- readRDS("ABCDv2.0.1_update_DEAP_baseline_sMRI.Rds")
pcs<-read.table(paste(geno_dir,"preImputation_QC/intermediate_files/sampleQC/PopStr/","ABCD_release_3.0_QCed_EUR_6sd.pop_strat_outliers.pca.evec",sep=""),stringsAsFactors = FALSE)
colnames(pcs)<-c("id",paste("PC",1:100,sep=""),"pheno")
pcs$subjectid<-strsplit(pcs$id,":") %>% sapply("[[",2) %>% factor()
#
d<-merge(dat,info,by.x="subjectid","abcd.id_redcap")
d<-merge(d,pcs,by="subjectid")

#---------------------------------------------------
# define variables and covariates
# covariates
pc_vars=paste("PC",1:10,sep="")
covs_def=c("sex","age",
           "Axiom_Plate",
           pc_vars
           )
vars<-c(covs_def,"abcd_site","mri_info_device.serial.number")
#
qcovs<-c("age",pc_vars)
ccovs<-vars[!vars %in% qcovs]

# save covariates to include in reml analysis
d1<-d[,c("RUID","subjectid",vars)]
d1<-d1 %>% rename(FID=RUID,IID=subjectid)

write.table(d1[,c("FID","IID",qcovs)],file=paste(base_dir,"covariates_q.txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE)
write.table(d1[,c("FID","IID",ccovs)],file=paste(base_dir,"covariates_c.txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE)
write.table(d1[,c("FID","IID",ccovs[c(1:2,3)])],file=paste(base_dir,"covariates_c_behaviour.txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE)
write.table(d1[,c("FID","IID",ccovs[c(1:2,4)])],file=paste(base_dir,"covariates_c_scanner.txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE)

