# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
i=args[2]
#
config_file="reading.config" # parameters for this run
i=1 # args[2]
#---------------------------------------------------
# define libraries
# define libraries
library(dplyr); library(tidyr)
library(data.table)
library(gamm4) #
library(lme4)
library(MuMIn)
library(psych)
library(lmerTest)

#---------------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
source(paste0(scripts_dir,"general_functions.R"))
#---------------------
# source config file to define parameters for this run
source(paste(scripts_dir,config_file,sep=""))
# source(config_file)

#---------------------
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",sep="")
input_dir=paste(pheno_dir,"baseline","_",atlas,"/",sep="") %>% gsub(" ","",.)
working_dir=paste(pheno_dir,atlas,sep="")
if(!dir.exists(working_dir)){dir.create(working_dir)}
setwd(working_dir)
#---------------------------------------------------
# model to run, as baseline
models_table<-read.csv(paste(input_dir,"trim",trim_val,"/tables/baseline_models_table.csv",sep=""),stringsAsFactors = FALSE)

m=models_table[i,"model_name"] %>% as.character()
model_b=models_table[i,"formula"] %>% as.character() # baseline model
# formulas for global measure adjusted models
model_b_tCSA<-models_table %>% filter(model_name=="m0c") %>% pull(formula)
model_b_mCT<-models_table %>% filter(model_name=="m0d") %>% pull(formula)

#---------------------------------------------------
# define output directories (dependent on cov adjustement)
if(length(grep("0",m))>0){
  out_dir=paste(working_dir,"/","baseline","/",sep="")
} else {
  cov1<-strsplit(model_b,"PC10") %>% sapply("[[",2) %>% strsplit(.," \\+ ") %>% sapply("[[",2)   
  out_dir=paste(working_dir,"/",cov1,"/",sep="")
  rm(cov1)
}

if (!dir.exists(out_dir)){dir.create(out_dir)}
if(!dir.exists(paste(out_dir,"/tables/",sep=""))){dir.create(paste(out_dir,"/tables/",sep=""))}
if(!dir.exists(paste(out_dir,"/figures/",sep=""))){dir.create(paste(out_dir,"/figures/",sep=""))}

#---------------------------------------------------
# read data
f=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.csv",sep="")
d<-read.csv(f)
d$rel_family_id <- as.factor(d$rel_family_id)
#---------------------------------------------------
# define rois
rois2follow <- read.table(paste(pheno_dir,atlas,"summary","tables","ROImeasures2follow.txt",sep="/"))

rois_rh<-gsub("_lh","_rh",rois2follow$V1) # only the ones that show assoc. with LH counterpart
rois_lh<- rois2follow$V1 %>% as.character()
rois<-c(rois_lh,rois_rh)
#---------------------------------------------------

#---------------------------------------------------
# Regression
# for each brain measure selected to follow up, run regression including covariates in baseline model, and get residuals
#---------------------------------------------------
d_res<-d %>% select(subjectid)

#---------------------------------------------------
## Run models

#-----------------------------------
for (bL in rois_lh) {
  bR<-gsub("lh","rh",bL)
  
  cat("Line:",bL,"\n")
  
  # define models, with brain measures as dv
  model_b0L<-gsub(dv,bL,model_b)
  model_b0R<-gsub("lh","rh",model_b0L)
  
  # global measure adjusted models, dependent on type of measure
  if(length(grep("area",bL))==1){
    model_b1L<-gsub(dv,bL,model_b_tCSA)
  } else if (length(grep("thick",bL))==1){
    model_b1L<-gsub(dv,bL,model_b_mCT)
  } else {
    model_b1L<-NULL
  }
  model_b1R<-gsub("lh","rh",model_b1L)
  
  
  cat("Baseline model (LH):",model_b0L,"\n")
  
  
  ## Run models
  # baseline model
  # m0.gam<-gamm4(as.formula(paste(dv,"~",paste(covs_def,collapse=" + "))),random=as.formula(paste("~",random)), data=d)
  m_b0L<-do.call("lmer",list (as.formula(model_b0L),data=d)) # LH
  m_b0R<-do.call("lmer",list (as.formula(model_b0R),data=d)) # RH
  # global measure adjusted models
  m_b1L<-do.call("lmer",list (as.formula(model_b1L),data=d)) # LH
  m_b1R<-do.call("lmer",list (as.formula(model_b1R),data=d)) # RH
  
  # anova(m_b0L,m_b1L)
  # anova(m_b0R,m_b1R)
  
  #---------------------------------------------------
  d_res[,paste0("residuals_",bL)]<-residuals(m_b0L)
  d_res[,paste0("residuals_",bR)]<-residuals(m_b0R)
  d_res[,paste0("residuals_globalAdj_",bL)]<-residuals(m_b1L)
  d_res[,paste0("residuals_globalAdj_",bR)]<-residuals(m_b1R)
  #---------------------------------------------------
  
  # clean intermediate
  rm(m_b0L,m_b0R,m_b1L,m_b1R)
  rm(model_b0L,model_b0R,model_b1L,model_b1R)
}


# save residuals
f2=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis_ROIs2follow_residuals.csv",sep="")
write.csv(d_res,file=f2,row.names=FALSE)

# repeat running gam and checking those residuals, just in case...


