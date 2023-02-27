# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
# config_file<-args[1]
# i=args[2]
#
config_file="reading.config" # parameters for this run
#i=2
#---------------------------------------------------
# define libraries
library(dplyr); library(tidyr)
library(data.table)
library(gamm4) #
library(lme4)
library(MuMIn)
library(psych)
library(lmerTest)
# to plot interactions
library(sjPlot)
library(sjmisc)
library(ggplot2)

# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
#---------------------
# source config file to define parameters for this run
source(paste(scripts_dir,config_file,sep=""))
# source(config_file)

#---------------------
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",sep="")
input_dir=paste(pheno_dir,"baseline","_",atlas,"/",sep="") %>% gsub(" ","",.)
working_dir=paste(pheno_dir,atlas,"/",sep="")
if(!dir.exists(working_dir)){dir.create(working_dir)}
setwd(working_dir)

#---------------------------------------------------
# read data
f=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.RDS",sep="")
d<-readRDS(f)
#---------------------------------------------------
# define rois
rois_lh<-read.table(paste(working_dir,"/summary/tables/ROImeasures2follow.txt",sep=""),header=FALSE) %>% pull(V1) %>% as.vector()

#---------------------------------------------------
models_table<-read.csv(paste(input_dir,"trim",trim_val,"/tables/baseline_models_table.csv",sep=""),stringsAsFactors = FALSE)

for (i in c(1,4:5)){
# model to run, as baseline
m=models_table[i,"model_name"] %>% as.character()
model_b=models_table[i,"formula"] %>% as.character() # baseline model

# define output directories (dependent on cov adjustement)
if(length(grep("0",m))>0){
  out_dir=paste(working_dir,"baseline","/",sep="")
} else {
  cov1<-strsplit(model_b,"PC10") %>% sapply("[[",2) %>% strsplit(.," \\+ ") %>% sapply("[[",2)   
  out_dir=paste(working_dir,cov1,"/",sep="")
  rm(cov1)
}

if (!dir.exists(out_dir)){dir.create(out_dir)}
if(!dir.exists(paste(out_dir,"/tables/",sep=""))){dir.create(paste(out_dir,"/tables/",sep=""))}
if(!dir.exists(paste(out_dir,"/figures/",sep=""))){dir.create(paste(out_dir,"/figures/",sep=""))}

#---------------------------------------------------
# Analysis
## re-run regressions including interactions between the regional measure and:
##    - sex

#---------------------------------------------------
## Run models

cat("Line:",i,"\n")
cat("Baseline model:",m,"\n")
cat(model_b,"\n")

# baseline model
m_b<-do.call("lmer",list (as.formula(model_b),data=d))

if(i!=1){
  # extract global term to evaluate interaction with
  global_var<-strsplit(model_b," \\+ ") %>% unlist() %>% tail(.,n=1)
  
  # m_b_effects_plot<-plot_model(m_b,show.values=TRUE) +
  #   theme_sjplot() +
  #   theme(legend.position="bottom")
  # ggsave(m_b_effects_plot,file=paste0(out_dir,"/figures/","interaction_modelEffects_",global_var,".png"))
  # 
  # m_b_pred_plot<-plot_model(m_b, type = "pred", terms = c(global_var)) +
  #   theme_sjplot()
  # ggsave(m_b_pred_plot,file=paste0(out_dir,"/figures/",global_var,".png"), height=5, width=5)
  
  
}

# RUN, one model per ROI

# define which rois to run analyses on
## if corrected for thickness -> only thickness
## if corrected for area -> only area
## if corrected for TBV or none  --> thickess and area
if (length(grep("thick",model_b))==1){
  rois2run<-rois_lh[grep("thick",rois_lh)]
} else if (length(grep("area",model_b)==1)){
  rois2run<-rois_lh[grep("area",rois_lh)]
} else {
  rois2run<-rois_lh
}
int_var="sex"
# only run if output does not exist
if (file.exists(paste(out_dir,"/tables/baseline_",m,"intSex_",atlas,"_","summary_t_table.csv",sep=""))==FALSE) # only run if output file does not exist 
{
  for (iv in rois2run){
    # if(iv!=global_var){
      # define test model
      model_t<-paste( model_b,"+ ",iv, sep=""    )
      m_t<-do.call("lmer",list (as.formula(model_t),data=d))
      # define model with intereaction
      model_int_t<-paste( model_b,
                      "+ ",
                      int_var, "*", iv, sep=""
      )
      m_int_t<-do.call("lmer",list (as.formula(model_int_t),data=d))
      # # check and save plots
      # m_effects_plot<-plot_model(m_t,show.values=TRUE) +
      #   theme_sjplot() +
      #   theme(legend.position="bottom")
      # m_int_effects_plot<-plot_model(m_int_t,show.values=TRUE) +
      #   theme_sjplot() +
      #   theme(legend.position="bottom")
      # m_int_plot<-plot_model(m_int_t, type = "pred", terms = c(int_var,iv)) +
      #   theme_sjplot()
      # ggsave(m_int_plot,file=paste0(out_dir,"/figures/","interaction_",int_var,"_",iv,".png"),width=5,height=5)
      # ggsave(m_int_effects_plot,file=paste0(out_dir,"/figures/","interaction_modelEffects_",int_var,"_",iv,".png"),width=5,height=5)
      # ggsave(m_effects_plot,file=paste0(out_dir,"/figures/","modelEffects_",int_var,"_",iv,".png"),width=5,height=5)
      # rm(m_int_plot,m_int_effects_plot,m_effects_plot); gc()
      #------------
      # summary table, parameter table - to get t-values
      sumt<-summary(m_int_t) %>% coef()
      int_term<-rownames(sumt)[NROW(sumt)]
      int_term2<- gsub("sexM","sex",int_term)
      sumt<-sumt[int_term,] %>% t %>% data.frame() %>% mutate(model=int_term)
      # ANOVA table, for all parameters
      anovat<-anova(m_int_t)[int_term2,] %>% data.frame() %>% mutate(model=int_term)
      # LRT  (REML)
      lrt<-anova(m_t,m_int_t) %>% data.frame()
      lrt$model<-c(paste("test",i,sep=""),int_term)
      # effect size, r2 for the fixed effects(r2m)
      effect_size<-data.frame(model=int_term,r2m=r.squaredGLMM(m_int_t)[1],delta_r2m=(r.squaredGLMM(m_int_t)[1]-r.squaredGLMM(m_t)[1]))
      #------------
      # save to objects across all rois
      if(exists("anova_table")){anova_table<-rbind(anova_table,anovat)} else {anova_table<-anovat}
      if(exists("sum_table")){sum_table<-rbind(sum_table,sumt)} else {sum_table<-sumt}
      if(exists("lrt_table")){lrt_table<-merge(lrt_table,lrt,all=TRUE,stringsAsFactors=FALSE)} else {lrt_table<-lrt }
      if(exists("effect_table")){effect_table<-merge(effect_table,effect_size,all=TRUE,stringsAsFactors=FALSE)} else {effect_table<-effect_size}
      
      # clean intermediate
      rm(m_t,model_t,
         model_int_t, m_int_t,
         sumt,anovat,lrt,effect_size)
    # }
  }
  # save summary tables
  write.csv(anova_table,file=paste(out_dir,"/tables/baseline_",m,"intSex_",atlas,"_","anova_table.csv",sep=""),row.names = FALSE)
  write.csv(sum_table,file=paste(out_dir,"/tables/baseline_",m,"intSex_",atlas,"_","summary_t_table.csv",sep=""),row.names = FALSE)
  write.csv(lrt_table,file=paste(out_dir,"/tables/baseline_",m,"intSex_",atlas,"_","LRT_table.csv",sep=""),row.names = FALSE)
  write.csv(effect_table,file=paste(out_dir,"/tables/baseline_",m,"intSex_",atlas,"_","effect_r2m_table.csv",sep=""),row.names = FALSE)
  # clean
  rm(anova_table,sum_table,lrt_table,effect_table)
}

rm(m_b)
#---------------------------------------------------
}