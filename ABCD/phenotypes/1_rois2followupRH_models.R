# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
i=args[2]
#
# config_file="reading.config" # parameters for this run
# i=1 # args[2]
#---------------------------------------------------
# define libraries
library(dplyr); library(tidyr)
library(data.table)
library(gamm4) #
library(lme4)
library(MuMIn)
library(psych)
library(lmerTest)

# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
#---------------------
# source config file to define parameters for this run
# source(paste(scripts_dir,config_file,sep=""))
source(config_file)

#---------------------
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
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
model_b<-gsub(".lh",".rh",model_b)

# define output directories (dependent on cov adjustement)
if(length(grep("0",m))>0){
  out_dir=paste(working_dir,"rois2followup_RH_","baseline","/",sep="")
} else {
  cov1<-strsplit(model_b,"PC10") %>% sapply("[[",2) %>% strsplit(.," \\+ ") %>% sapply("[[",2)   
  out_dir=paste(working_dir,cov1,"/",sep="")
  rm(cov1)
}


if (!dir.exists(out_dir)){dir.create(out_dir)}
if(!dir.exists(paste(out_dir,"/tables/",sep=""))){dir.create(paste(out_dir,"/tables/",sep=""))}
if(!dir.exists(paste(out_dir,"/figures/",sep=""))){dir.create(paste(out_dir,"/figures/",sep=""))}

#---------------------------------------------------
# read data
f=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.csv",sep="")
d<-read.csv(f)
d$rel_family_id<-as.factor(d$rel_family_id)
#---------------------------------------------------
# define rois
rois2follow <- read.table(paste(pheno_dir,atlas,"summary","tables","ROImeasures2follow.txt",sep="/"))
rois_rh<-gsub(".lh",".rh",rois2follow$V1) # only the ones that show assoc. with LH counterpart
#---------------------------------------------------

#---------------------------------------------------
## Run models

cat("Line:",i,"\n")
cat("Baseline model:",m,"\n")
cat(model_b,"\n")

# baseline model
m_b<-do.call("lmer",list (as.formula(model_b),data=d))
# RUN, one model per ROI

# define which rois to run analyses on
## if corrected for thickness -> only thickness
## if corrected for area -> only area
## if corrected for TBV or none  --> thickess and area
if (length(grep("thick",model_b))==1){
  rois2run<-rois_rh[grep("thick",rois_rh)]
} else if (length(grep("area",model_b)==1)){
  rois2run<-rois_rh[grep("area",rois_rh)]
} else {
  rois2run<-rois_rh
}

# only run if output does not exist
if (file.exists(paste(out_dir,"/tables/baseline_",m,"_",atlas,"_","summary_t_table.csv",sep=""))==FALSE) # only run if output file does not exist 
  {
  counter=1
  for (iv in rois2run){
    # define test model
    model_t<-paste( model_b, 
                    iv, sep=" + "
    )
    m_t<-do.call("lmer",list (as.formula(model_t),data=d))
    
    # summary table, parameter table - to get t-values
    sumt<-summary(m_t) %>% coef()
    sumt<-sumt[iv,] %>% t %>% data.frame() %>% mutate(model=iv)
    # ANOVA table, for all parameters
    anovat<-anova(m_t)[iv,] %>% data.frame() %>% mutate(model=iv)
    # LRT  (REML)
    lrt<-anova(m_b,m_t) %>% data.frame()
    lrt$model<-c(paste("baseline",i,sep=""),iv)
    # effect size, r2 for the fixed effects(r2m)
    effect_size<-data.frame(model=iv,r2m=r.squaredGLMM(m_t)[1],delta_r2m=(r.squaredGLMM(m_t)[1]-r.squaredGLMM(m_b)[1]))
    
    # save to objects across all rois
    if(exists("anova_table")){anova_table<-rbind(anova_table,anovat)} else {anova_table<-anovat}
    if(exists("sum_table")){sum_table<-rbind(sum_table,sumt)} else {sum_table<-sumt}
    if(exists("lrt_table")){lrt_table<-merge(lrt_table,lrt,all=TRUE,stringsAsFactors=FALSE)} else {lrt_table<-lrt }
    if(exists("effect_table")){effect_table<-merge(effect_table,effect_size,all=TRUE,stringsAsFactors=FALSE)} else {effect_table<-effect_size}
    
    # clean intermediate
    rm(m_t,model_t,m_t,sumt,anovat,lrt,effect_size)
    
    # echo some runs
    counter=counter+1
    if(counter %in% seq(0,length(rois_rh),by=25)){
      cat("Running phenotype number:", counter,"\n")
      cat("Phenotype: ", iv,"\n")
      cat("--------------------------------------------------------------")
    }
  }
  rm(counter)
  # save summary tables
  write.csv(anova_table,file=paste(out_dir,"/tables/baseline_",m,"_",atlas,"_","anova_table.csv",sep=""),row.names = FALSE)
  write.csv(sum_table,file=paste(out_dir,"/tables/baseline_",m,"_",atlas,"_","summary_t_table.csv",sep=""),row.names = FALSE)
  write.csv(lrt_table,file=paste(out_dir,"/tables/baseline_",m,"_",atlas,"_","LRT_table.csv",sep=""),row.names = FALSE)
  write.csv(effect_table,file=paste(out_dir,"/tables/baseline_",m,"_",atlas,"_","effect_r2m_table.csv",sep=""),row.names = FALSE)
  # clean
  rm(anova_table,sum_table,lrt_table,effect_table)

}
#---------------------------------------------------