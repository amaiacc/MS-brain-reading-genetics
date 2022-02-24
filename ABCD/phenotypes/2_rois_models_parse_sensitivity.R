# clean workspace
rm(list=ls())
library(ggseg)
library(ggsegDesterieux)
#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/reading.config" # parameters for this run
#---------------------------------------------------

# define libraries
library(dplyr); library(tidyr)
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
# source config file to define parameters for this run
source(config_file)
#---------------------
if(atlas=="destrieux"){
  target_labels<-data.frame(target=ggseg(atlas=desterieux)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
} else if (atlas=="desikan"){
  target_labels<-data.frame(target=ggseg(atlas=dkt)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
}

base_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",atlas,"/",sep="")
setwd(base_dir)

#---------------------
hemi="lh"
# save tables
working_dir=paste(base_dir,"summary","/",sep="")
out_dir=paste(working_dir,"/tables/",sep="")
if(!dir.exists(working_dir)){dir.create(working_dir)}
if(!dir.exists(out_dir)){dir.create(out_dir)}

rois2follow<-read.table(paste(out_dir,"ROImeasures2follow.txt",sep=""),header=FALSE) %>% pull(V1) %>% as.vector()


## sensitivity analyses controlling for cognitive covariates
for(cov1 in c(cov11,cov12)){
  
  ttables<-c( paste(paste(cov1,"/tables/",sep=""),list.files(paste(cov1,"/tables/",sep=""),pattern="t_table"),sep="/"))
  ttables<-ttables[grep("b_|b2.3_",ttables,invert=TRUE)]

  for (f in ttables){
    i=strsplit(f,"_m") %>% sapply("[[",2) %>% strsplit(.,"_") %>% sapply("[[",1)
    
    # build c (model name), sequentially
    c<-NULL
    if ( length(grep("^1",i))==1 ){
      c<-paste(c,cov1,sep=" + ")
    }
    if ( length(grep("^2",i))==1 ){
      c<-paste(c,cov1,cov2,sep=" + ")
    }
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
  
    # read data
    tvals<-read.csv(f)
    lrt<-read.csv(gsub("_summary_t_","_LRT_",f))
    t<-merge(tvals,lrt,by="model",suffixes=c(".t",".LRT"))
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
    t$Pfdr<-p.adjust(t$P,method = "fdr")
    t$Pbonf<-p.adjust(t$P,method = "bonferroni")
    t$adj<-c
    t$f<-f
    # combine data
    if(exists("d")){d<-rbind(d,t)} else {d<-t}
    # clean intermediate
    rm(c,t,i)
    rm(tvals,lrt)
  }
  rm(f)
}
rm(cov1)


# keep track of two types of models
d <- d %>% mutate(regModel=if_else(grepl("Mean|Total",adj),"Model2","Model1"))

#
# which regions are associated with reading? fdr corrected
d<-d %>% select(regModel,adj,model,measure,region,hemi,label,
                Estimate,StdError,df,t.value,P,Pfdr,Pbonf,
                AIC,BIC,logLik,deviance,Chisq,Df,Pr..Chisq.) %>%
  rename(indep.variable=model)

d$adj <- d$adj %>% gsub("TotalAreaLH","totalCSAlh",.) %>% 
  gsub("MeanThicknessLH","meanCTlh",.) %>%
  gsub("nihtbx_|pea_","",.) %>% gsub("_uncorrected|_tss","",.)
# subset only to the rois to follow up, to check the rest of sensitivity analyses and genetic analyses
d_subset<-subset(d,indep.variable %in% rois2follow) %>% arrange(indep.variable,adj) %>% filter(adj!="TBV"&adj!="TBV^2/3")


# save tables
write.csv(d,paste(out_dir,"sensitivity_all_",atlas,"_summary_t_LRT_table.csv",sep=""),row.names=FALSE)
write.csv(d_subset,paste(out_dir,"sensitivity_sig_",atlas,"_summary_t_LRT_table.csv",sep=""),row.names=FALSE)


