# perform mediation analyses, to explore the relationship between:
## dv: reading
## iv: PGS (CP)
## mediation: CSA (CSAtotal and CSTstg)
# is the association between reading and PGS-CP mediated through CSA?
#---------------------

# clean workspace
rm(list=ls())

#---------------------
# args <- commandArgs(TRUE)
# config_file<-args[1]
# i=args[2]
#
config_file="reading.config" # parameters for this run
#---------------------------------------------------
# define libraries
library(dplyr); library(tidyr)
library(data.table)
# library(gamm4) #
library(lme4)
# library(MuMIn)
library(psych)
# library(lmerTest)
library(mediation)
#---------------------------------------------------

#---------------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {
  dir="F:/projects/"
  } else {
    if(length(grep("cajal",Sys.info()['nodename']))==1){
      dir="/bcbl/home/home_a-f/acarrion/projects/"
    } else {
      dir="/export/home/acarrion/acarrion/projects/"
    }
  }
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
#---------------------
# source to get custom functions for mediation
source(paste0(scripts_dir,"functions_mediation.R"))

# source config file to define parameters for this run
source(paste(scripts_dir,config_file,sep=""))
# source(config_file)
i=1 # args[2]
s="eur6sd_unrelated" # s="all_unrelated"

#---------------------
pgs_dir=paste(dir,"/resources/datasets/ABCD/data/working_data/genotyping/prsice/cognitive/QC/",sep="")
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",sep="")
input_dir=paste(pheno_dir,"baseline","_",atlas,"/",sep="") %>% gsub(" ","",.)
working_dir=paste(pheno_dir,atlas,"/",sep="")
if(!dir.exists(working_dir)){dir.create(working_dir)}
setwd(working_dir)
#---------------------------------------------------
# define output directories
out_dir=paste(working_dir,"mediation","/",sep="")
if (!dir.exists(out_dir)){dir.create(out_dir)}
if(!dir.exists(paste(out_dir,"/tables/",sep=""))){dir.create(paste(out_dir,"/tables/",sep=""))}
if(!dir.exists(paste(out_dir,"/figures/",sep=""))){dir.create(paste(out_dir,"/figures/",sep=""))}
#---------------------------------------------------
# read data
f=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.csv",sep="")
d<-read.csv(f)
d$rel_family_id<-as.factor(d$rel_family_id)
#---------------------------------------------------
# read polygenic score for PC (thr=0.05)
pgs<-read.csv(paste0(pgs_dir,"PGS_PCs_cognitive.csv"))[,c("id","CP_0.05")]
# combine
d<-merge(d,pgs,by.x="subjectid",by.y="id")
rm(pgs)
#---------------------------------------------------
# z-scores for all the quantiative dependent and independent variables
d<-d %>% mutate_if(is.numeric,scale)
d$rel_family_id<-factor(d$rel_family_id)
# make sure class is numeric (not matrix as output from the scaling) otherwise mediation complains
d[,dv]<-as.numeric(d[,dv])
#---------------------------------------------------
# define rois
rois<-colnames(d)[grep(atlas,colnames(d))]
rois<-rois[grep("area",rois)]
rois_lh<-rois[grep("lh",rois)]
controls<-c("smri_thick_cort.destrieux_mean.lh",
            "smri_thick_cort.destrieux_g.and.s.occipital.inf.lh",
            "smri_thick_cort.destrieux_g.oc.temp.lat.fusifor.lh",
            "smri_area_cort.destrieux_g.parietal.sup.lh",
            "smri_area_cort.destrieux_g.temp.sup.lateral.rh",
            "smri_area_cort.destrieux_total.rh") # add as neg control...
#
cogn_covs0=c("nihtbx_fluidcomp_uncorrected","pea_wiscv_tss")
cogn_covs1=c(cogn_covs0,"nihtbx_picvocab_uncorrected")
#---------------------------------------------------
# Define models to run
models_table<-read.csv(paste(input_dir,"trim",trim_val,"/tables/baseline_models_table.csv",sep=""),stringsAsFactors = FALSE)
m=models_table[i,"model_name"] %>% as.character()

# baseline model
model_b=models_table[i,"formula"] %>% as.character() 

# get part of model adjusting for random structure and covariates
model_covs<-gsub(dv,"",model_b)
# get covs as list
covs=strsplit(model_covs,") \\+ ") %>% sapply("[[",2) %>% strsplit(.," \\+ ") %>% unlist()
random=strsplit(model_covs," \\+ ") %>% sapply("[[",1) %>% unlist() %>% gsub("~ ","",.)
#---------------------------------------------------
# MEDIATION
# Mediation of CSA on CP effect on reading
## dv already defined by the config file, it's constant
#---------------------------------------------------
CSAt<-rois_lh[grep("total",rois_lh)]
CSAstg<-rois_lh[grep("g.temp.sup.lateral.lh",rois_lh)]
#---------------------------------------------------
# dv: reading
# iv: pgs
iv<- "CP_0.05"
#---------------------------------------------------
# mediator: CSA (total)
mediator=CSAt
#
read_med_CSAt<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=covs,random=random,d=d)
read_med_CSAt_adjCSAstg<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,CSAstg),random=random,d=d)
read_med_CSAt_adjCogn0<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,cogn_covs0),random=random,d=d)
read_med_CSAt_adjCogn1<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,cogn_covs1),random=random,d=d)

rm(mediator)
#---------------------------------------------------
# mediator: CSA (stg)
mediator=CSAstg
#
read_med_CSAstg<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=covs,random=random,d=d)
read_med_CSAstg_adjCSAt<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,CSAt),random=random,d=d)
read_med_CSAstg_adjCogn0<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,cogn_covs0),random=random,d=d)
read_med_CSAstg_adjCogn1<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,cogn_covs1),random=random,d=d)

rm(mediator)

#---------------------------------------------------
# summarize all runs for reading
read_mediation_all<-lapply(ls(pattern="read_med_"),function(x){
  get(x) %>% get_mediation_table() %>% 
    mutate(mediation_name=x)
}) %>% do.call("rbind",.) %>%
  mutate(adjustement=gsub(paste(covs,collapse=","),"",covariates) %>% gsub("^,","",.)) 

#---------------------------------------------------
# clean tables to include only relevant columns...
read_clean <- read_mediation_all %>% distinct() %>%
  dplyr::select(dv,iv,mediator,adjustement,stat,Estimate,CIupper,CIlower,p,adjustement,mediation_name)
#---------------------------------------------------

#---------------------------------------------------
# save results from the mediation analyses
write.csv(read_clean,paste0(out_dir,"/tables/reading_PGS_mediation_v2_summary_MC10000.csv"),row.names=FALSE)

save(list = ls(all.names = TRUE,pattern="read_med"), file = paste0(out_dir,"reading_PGS_mediation_v2_MC10000.RData"))
#---------------------------------------------------