# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/reading.config" # parameters for this run
#---------------------------------------------------

# define libraries
library(dplyr); library(tidyr)
library(ggseg)
library(ggsegDefaultExtra)
library(ggsegDesterieux)
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())

# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
# source config file to define parameters for this run
# source(paste(scripts_dir,config_file,sep="")) # if hard-coded
source(config_file) # as parameter
#---------------------
if(atlas=="destrieux"){
  target_labels<-data.frame(target=ggseg(atlas=desterieux)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
} else if (atlas=="desikan"){
  target_labels<-data.frame(target=ggseg(atlas=dkt)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
}

base_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",atlas,"/",sep="")
setwd(base_dir)

# save tables
working_dir=paste(base_dir,"summary","/",sep="")
out_dir=paste(working_dir,"/tables/",sep="")
if(!dir.exists(working_dir)){dir.create(working_dir)}
if(!dir.exists(out_dir)){dir.create(out_dir)}

#---------------------
hemi="lh"

# baseline models
t_tables<-c( paste(paste("baseline","/tables/",sep=""),
                  list.files(paste("baseline","/tables/",sep=""),pattern="t_table"),sep="/"))
t_tables<-t_tables[grep("intSex|intAge",t_tables)]


for (f in t_tables){
  i=strsplit(f,"_m") %>% sapply("[[",2) %>% strsplit(.,"_") %>% sapply("[[",1)
  # build c (model name), sequentially
  c<-NULL
  if ( length(grep("b$",i))==1 ){
    c<-paste(c,"TBV",sep=" + ")
  }
  if ( length(grep("2.3$",i))==1 ){
    c<-paste(c,"TBV^2/3",sep=" + ")
  }
  if ( length(grep("c$",i))==1 ){
    c<-paste(c,"TotalAreaLH",sep=" + ") 
  }
  if ( length(grep("d$",i))==1 ){
    c<-paste(c,"MeanThicknessLH",sep=" + ")
  }
  if(is.null(c)){c<-""}
  c<-gsub(" \\+ ","",c)
  # read data
  tvals<-read.csv(f)
  lrt<-read.csv(gsub("_summary_t_","_LRT_",f))
  t<-merge(tvals,lrt,by="model",suffixes=c(".t",".LRT"))
  t$file<-f
  # create additional columns for organization/visualization
  t$hemi<-NA;
  t$hemi[grep(".lh",t$model)]<-"lh"
  t$hemi[grep(".rh",t$model)]<-"rh"
  t$measure<-strsplit(as.character(t$model),"_") %>% sapply("[[",2)
  t$region<-strsplit(as.character(t$model),atlas) %>% sapply("[[",2) %>% gsub("^_","",.) %>% gsub("\\.rh|\\.lh","",.) 
  t$matchLabel<-paste(t$hemi,t$region,sep="_")  %>% gsub("\\.","_",.)
  t<-merge(t,target_labels,by="matchLabel",all.x=TRUE) %>% mutate(label=as.factor(target)) %>% select(-matchLabel,-target) 
  t<-t %>% rename(P= Pr...t..,StdError=Std..Error) %>% mutate(minuslog10P=(-log10(P)))
  t$label<-as.factor(t$label)
  # t<- t %>% filter(region!="total" & region!="mean")
  t$adj<-c
  t$Pfdr<-p.adjust(t$Pr..Chisq.,method = "fdr")
  t$Pbonf<-p.adjust(t$Pr..Chisq.,method = "bonferroni")
  # combine data
  if(exists("d")){d<-rbind(d,t)} else {d<-t}
  # clean intermediate
  rm(c,t,i)
  rm(tvals,lrt)
}
rm(f)

# keep track of two types of models
d <- d %>% mutate(regModel=if_else(adj=="","Model1","Model2"))
# remove model2 with global measures, which is the same as prev.
d <- d %>% filter(!(regModel=="Model2" & (region=="total"|region=="mean")))


# which regions are associated with reading? fdr corrected - for all tests included here...
d<-d %>% select(file,regModel,adj,model,measure,region,hemi,label,
                Estimate,StdError,df,t.value,P,Pfdr,Pbonf,
                AIC,BIC,logLik,deviance,Chisq,Df,Pr..Chisq.) %>%
  rename(indep.variable=model)

d$adj[d$adj==""]<-"none"
d$adj<-gsub(" \\+ ","",d$adj)
d$adj %>% table()
d$adj<-factor(d$adj,c("none","TotalAreaLH","MeanThicknessLH","TBV","TBV^2/3"))


# save files
write.csv(d,paste(out_dir,"baseline_rois2followup_sigAgeSex_",atlas,"_summary_t_LRT_table.csv",sep=""),row.names=FALSE)

