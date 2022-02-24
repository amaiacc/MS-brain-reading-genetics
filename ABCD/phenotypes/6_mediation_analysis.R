# perform mediation analyses, to explore the relationship between:
## dv: reading
## iv: CSAtemp
## mediation: CSAtotal
# is the association between reading and temporal CSA mediated through total CSA?
# or the other way around!?
# step2: include PGS into models...
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
library(gamm4) #
library(lme4)
library(MuMIn)
library(psych)
# library(lmerTest)
library(mediation)
#---------------------------------------------------

#---------------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
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
cogn_covs=c("nihtbx_fluidcomp_uncorrected","pea_wiscv_tss")
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
# Mediation of CP on CSA effect on reading
## dv already defined by the config file, it's constant
#---------------------------------------------------
CSAt<-rois_lh[grep("total",rois_lh)]
CSAstg<-rois_lh[grep("g.temp.sup.lateral.lh",rois_lh)]
#---------------------------------------------------
# dv: reading
# iv: pgs
mediator<- "CP_0.05"
#---------------------------------------------------
# iv: CSA (total)
iv=CSAt
#
read_med_CSAt<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=covs,random=random,d=d)
read_med_CSAt_adjCogn<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,cogn_covs),random=random,d=d)
read_med_CSAt_adjCSAstg<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,CSAstg),random=random,d=d)
read_med_CSAt_adjCSAtRH<-run_mediation(dv=dv,iv=iv,mediator=mediator,
                                         covs=c(covs,"smri_area_cort.destrieux_total.rh"),
                                         random=random,d=d)
rm(iv)
#---------------------------------------------------
# iv: CSA (stg)
iv=CSAstg
#
read_med_CSAstg<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=covs,random=random,d=d)
read_med_CSAstg_adjCogn<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,cogn_covs),random=random,d=d)
read_med_CSAstg_adjCSAt<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=c(covs,CSAt),random=random,d=d)
read_med_CSAstg_adjCSAstgRH<-run_mediation(dv=dv,iv=iv,mediator=mediator,
                                       covs=c(covs,"smri_area_cort.destrieux_g.temp.sup.lateral.rh"),
                                       random=random,d=d)

rm(iv)

#---------------------------------------------------
# mediation RH

# total
iv=gsub("lh","rh",CSAt)
read_med_CSAtRH<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=covs,random=random,d=d)
rm(iv)

# stg
iv=gsub("lh","rh",CSAstg)
read_med_CSAstgRH<-run_mediation(dv=dv,iv=iv,mediator=mediator,covs=covs,random=random,d=d)
read_med_CSAstgRH_adjCSAstg<-run_mediation(dv=dv,iv=iv,mediator=mediator,
                                 covs=c(covs,"smri_area_cort.destrieux_g.temp.sup.lateral.lh"),
                                 random=random,d=d)
rm(iv)

#---------------------------------------------------
# DV: fluid IQ measures
# iv: total
iv=CSAt
fiq1_med_CSAt<-run_mediation(dv=cogn_covs[1],iv=iv,mediator=mediator,covs=covs,random=random,d=d)
fiq1_med_CSAt_adjRead<-run_mediation(dv=cogn_covs[1],iv=iv,mediator=mediator,covs=c(covs,dv),random=random,d=d)
fiq2_med_CSAt<-run_mediation(dv=cogn_covs[2],iv=iv,mediator=mediator,covs=covs,random=random,d=d)
fiq2_med_CSAt_adjRead<-run_mediation(dv=cogn_covs[2],iv=iv,mediator=mediator,covs=c(covs,dv),random=random,d=d)
rm(iv)
# iv: stg
iv=CSAstg
fiq1_med_CSAstg<-run_mediation(dv=cogn_covs[1],iv=iv,mediator=mediator,covs=covs,random=random,d=d)
fiq1_med_CSAstg_adjRead<-run_mediation(dv=cogn_covs[1],iv=iv,mediator=mediator,covs=c(covs,dv),random=random,d=d)
fiq2_med_CSAstg<-run_mediation(dv=cogn_covs[2],iv=iv,mediator=mediator,covs=covs,random=random,d=d)
fiq2_med_CSAstg_adjRead<-run_mediation(dv=cogn_covs[2],iv=iv,mediator=mediator,covs=c(covs,dv),random=random,d=d)
rm(iv)
#---------------------------------------------------
# summarize all runs for reading
read_mediation_all<-lapply(ls(pattern="read_med_"),function(x){
  get(x) %>% get_mediation_table() %>% 
    mutate(mediation_name=x)
}) %>% do.call("rbind",.) %>%
  mutate(adjustement=gsub(paste(covs,collapse=","),"",covariates) %>% gsub("^,","",.)) 

# iq
fiq_mediation_all<-lapply(ls(pattern="fiq1_med_|fiq2_med_"),function(x){
  get(x) %>% get_mediation_table() %>% 
    mutate(mediation_name=x)
}) %>% do.call("rbind",.)  %>%
  mutate(adjustement=gsub(paste(covs,collapse=","),"",covariates) %>% gsub("^,","",.))


#---------------------------------------------------
# clean tables to include only relevant columns...
read_clean <- read_mediation_all %>% distinct() %>%
  dplyr::select(dv,iv,mediator,adjustement,stat,Estimate,CIupper,CIlower,p,adjustement,mediation_name)
iq_clean <- fiq_mediation_all %>% distinct() %>%
  dplyr::select(dv,iv,mediator,adjustement,stat,Estimate,CIupper,CIlower,p,adjustement,mediation_name)
#---------------------------------------------------

#---------------------------------------------------
# save results from the mediation analyses
write.csv(read_clean,paste0(out_dir,"/tables/reading_PGS_mediation2_summary.csv"),row.names=FALSE)
write.csv(iq_clean,paste0(out_dir,"/tables/fiq_PGS_mediation2_summary.csv"),row.names=FALSE)
save(list = ls(all.names = TRUE,pattern="read_med"), file = paste0(out_dir,"reading_PGS_mediation2.RData"))
save(list = ls(all.names = TRUE,pattern="fiq1_med|fiq2_med"), file = paste0(out_dir,"fiq_PGS_mediation2.RData"))
#---------------------------------------------------