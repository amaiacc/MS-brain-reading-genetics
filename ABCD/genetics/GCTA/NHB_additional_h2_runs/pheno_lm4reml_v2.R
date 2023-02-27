# clean workspace
rm(list=ls())

# parameters
args <- commandArgs(TRUE)
args<-c("ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned_rm05_adj",
        "smri_area_cort.destrieux_total",
        "ABCDv3_DEAP"
)

root<-args[1]
pheno<-args[2]
pheno_v<-args[3]

# define libraries
library(dplyr); library(tidyr)
library(lme4)

# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
pgeno_dir=paste(dir,"resources/datasets/ABCD/data/primary_data/ABCDgenomicsDataV30/genomics_sample03/ABCD_genotype/",sep="")
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/",sep="")

pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",sep="")
setwd(pheno_dir)
base_dir=paste(geno_dir,"/GCTA/pheno_v2/",sep="")
if (!dir.exists(base_dir)){dir.create(base_dir)}

# read data
ids<-read.table(paste(geno_dir,"/GCTA/grm/",root,".grm.id",sep=""))
# baseline info/phenotypes
s="eur6sd_unrelated"
atlas="destrieux"
trim_val="0.02"
input_dir=paste(pheno_dir,"/",pheno_v,"/",s,"/","reading","/","baseline","_",atlas,"/","trim",trim_val,"/",sep="") %>% gsub(" ","",.)
f=paste(input_dir,pheno_v,"_baseline_sMRI","_","reading","_","trim",trim_val,"_","scaled4analysis.csv",sep="")
d<-read.csv(f)
d$rel_family_id<-as.factor(d$rel_family_id)
#
d<-merge(d,ids,by.x="subjectid",by.y="V2",all=TRUE)
#
if(class(d[,pheno])!="numeric"){
  d[,pheno]<-as.character(d[,pheno]) %>% as.numeric()
}

# exclude Axiom Plate 416
## In our QC process, we found the plate 461 to be especially problematic. 
## We recommend not using genotype data from plate 461, or at least including plate number as a covariate.
# d<-subset(d,Axiom_Plate!="416")
#---------------------------------------------------
# define variables and covariates
# covariates
pc_vars=paste("PC",1:10,sep="")
covs_def=c("sex","age",
           # "Axiom_Plate",
           pc_vars
)

#
adjcols <- colnames(d)[grep("smri_thick_cort.destrieux_mean|smri_area_cort.destrieux_total",colnames(d))]
# select specific variable to adjust for, if any
if ( length(grep("smri",pheno))>0 & length(grep("total",pheno))==0 ) {
  m = strsplit(pheno,"_") %>% sapply("[[",2)
  if (length(grep("rh",pheno))>0){
    hemi="rh"
  } else if (length(grep("lh",pheno))>0){
    hemi="lh"
  }
  adjcov=adjcols[grep(paste(m,hemi,sep=".*."),adjcols)]
} else {
  adjcov=NULL
}

## all variables to consider
vars<-c(covs_def,pheno,"abcd_site","mri_info_device.serial.number")
#---------------------------------------------------
# define complete data
d<-d[complete.cases(d[,c(vars,adjcov)]),]

print(adjcov)

# scale numeric variables
d<-d %>% mutate_if(is.numeric,scale)

# define random part
if((grep("smri",pheno)%>% length > 0)==TRUE){
  random<-"(1|mri_info_device.serial.number)"
} else {
  random<-"(1|abcd_site)"
  
}

# run model
model_b0<-paste(pheno,"~",paste(c(random,covs_def),collapse=" + ",sep=" "),sep=" ")
print(model_b0)
m0<-do.call("lmer",list (as.formula(model_b0),data=d))

print("so far so good\n")
# save residuals
d2<-d[,c("V1","subjectid",pheno)]
d2[,c(paste(pheno,"_residuals",sep=""))]<-residuals(m0)
d2<-d2 %>% rename(FID=V1,IID=subjectid)

if (!file.exists(paste(base_dir,"residuals_",pheno,".txt",sep=""))){
 write.table(d2,file=paste(base_dir,"residuals_",pheno,".txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE)
}

# check phenotype, and run adjusted model if necessary 
## (i.e. imaging derived phenotype, not global)
if (!is.null(adjcov)){
  model_b1<-paste(pheno,"~",paste(c(random,covs_def,adjcov),collapse=" + ",sep=" "),sep=" ")
  print(model_b1)
  m1<-do.call("lmer",list (as.formula(model_b1),data=d))
  
print("so far so good\n")
# save residuals
d3<-d[,c("V1","subjectid",pheno)]
d3[,c(paste(pheno,"_residuals_adjGlobal",sep=""))]<-residuals(m1)
d3<-d3 %>% rename(FID=V1,IID=subjectid)

if (!file.exists(paste(base_dir,"residuals_",paste0("adjGlobal_",m,".",hemi),"_",pheno,".txt",sep=""))){
  write.table(d3,file=paste(base_dir,"residuals_",paste0("adjGlobal_",m,".",hemi),"_",pheno,".txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE)
}
}



