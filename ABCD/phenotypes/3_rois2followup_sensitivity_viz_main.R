# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/reading.config" # parameters for this run
# config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/picvocab.config" # parameters for this run
hemi="rh"
#---------------------------------------------------

# custom function to build c (model name), sequentially
model_name_build<-function(i,cov1,cov2){
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
  return(c)
}

#---------------------------------------------------
# define libraries
library(dplyr); library(tidyr)
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
#
library(ggseg)
# library(ggsegExtra)
library(ggsegDesterieux)
#library(ggsegDefaultExtra)
# library(ggseg3d)


# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
# source config file to define parameters for this run
source(config_file)
#---------------------
base_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",atlas,"/",sep="")
setwd(base_dir)
#---------------------
if(atlas=="destrieux"){
  target_labels<-data.frame(target=desterieux$data$label,
                            region=desterieux$data$region,
                            stringsAsFactors = FALSE) %>% 
    mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.)) %>% distinct() 
} else if (atlas=="desikan"){
  target_labels<-data.frame(target=dk$data$label,
                            region=dk$data$region,
                            stringsAsFactors = FALSE) %>% 
    mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.)) %>% distinct()
}


# get ROIs2follow, for main figures
rois2follow <-read.table("summary/tables/ROImeasures2follow.txt") %>% rename(model=V1) %>% select(model)


## sensitivity analyses controlling for cognitive covariates
for(cov1 in c(cov11,cov12)){
  
  ttables<-c( paste(paste("baseline","/tables/",sep=""),list.files(paste("baseline","/tables/",sep=""),pattern="t_table"),sep="/"),
              paste(paste(cov1,"/tables/",sep=""),list.files(paste(cov1,"/tables/",sep=""),pattern="t_table"),sep="/")
  )
  
  ttables<-ttables[grep("b_|b2.3_",ttables,invert=TRUE)] # do not include TBV adjusted - very similar to CSA and CT adjusted
  
  for (f in ttables){
    
    #
    i=strsplit(f,"_m") %>% sapply("[[",2) %>% strsplit(.,"_") %>% sapply("[[",1)
    
    # build c (model name), sequentially
    c<-model_name_build(i,cov1=cov1,cov2=cov2)
    
    # read data
    t<-read.csv(f)
    # create additional columns for organization/visualization
    t$hemi<-NA;
    t$hemi[grep(".lh",t$model)]<-"lh"
    t$hemi[grep(".rh",t$model)]<-"rh"
    t$measure<-strsplit(as.character(t$model),"_") %>% sapply("[[",2)
    t$region<-strsplit(as.character(t$model),atlas) %>% sapply("[[",2) %>% gsub("^_","",.) %>% gsub("\\.rh|\\.lh","",.) 
    t$matchLabel<-paste(t$hemi,t$region,sep="_")  %>% gsub("\\.","_",.)
    t<-merge(t,target_labels,by="matchLabel",all.x=TRUE,suffixes=c(".name","")) %>% 
      mutate(label=as.factor(target)) %>% select(-matchLabel,-target) 
    t<-t %>% rename(P= Pr...t..,StdError=Std..Error) %>% mutate(minuslog10P=(-log10(P)))
    t$label<-as.factor(t$label)
    # t<- t %>% filter(region!="total" & region!="mean")
    t$Pfdr<-p.adjust(t$P,method = "fdr")
    t$Pbonf<-p.adjust(t$P,method = "bonferroni")
    t$adj<-c
    # combine data
    if(exists("d")){d<-rbind(d,t)} else {d<-t}
    #clean intermediate
    rm(c,i)

  }
  rm(t)
}
rm(cov1)
#
d$region[is.na(d$region)]<-d$region.name[is.na(d$region)]

# select only rois2follow for plots
d$rois2follow<-FALSE
d$rois2follow[d$model %in% rois2follow$model ]<-TRUE

# new names for adjustement: i.e. adjGlobal or not
d$adjGlobal<-""
d$adjGlobal[grep("TotalAreaLH|MeanThicknessLH",d$adj)]<-"adjGlobal"
d$adjCognitive<-gsub(" \\+ TotalAreaLH| \\+ MeanThicknessLH","",d$adj) %>%
  gsub("nihtbx_|_uncorrected|pea_|_tss","",.)
#---------------------

#---------------------
# visualize estimates, across runs
library(RColorBrewer)

# get order of estimates, from the non adjusted LM:
betas_order_csa <- d %>% filter(adj==""&measure=="area") %>% arrange(Estimate)
betas_order_ct <- d %>% filter(adj==""&measure=="thick") %>% arrange(Estimate)

d_csa <- d %>% filter(rois2follow==TRUE) %>% filter(measure=="area") %>% filter(!(region=="total"&adj!="")) %>% 
  mutate(region=factor(region,levels=betas_order_csa %>% pull(region) %>% unique())) %>% 
  mutate(model=factor(model,levels=betas_order_csa %>% pull(model) %>% unique()))

d_ct<-d %>% filter(rois2follow==TRUE) %>% filter(measure=="thick")%>% filter(!(region=="mean"&adj!="")) %>% 
  mutate(region=factor(region,levels=betas_order_ct %>% pull(region) %>% unique())) %>% 
  mutate(model=factor(model,levels=betas_order_ct %>% pull(model) %>% unique()))

csa_betas <- ggplot(data=d_csa,
                    aes(x=Estimate,y=region,color=adjGlobal,shape=adjGlobal)) + 
  geom_vline(xintercept = 0,color="black",alpha=0.7,linetype="solid") +
  theme(axis.text.y=
          element_text(face=if_else(levels(d_csa$model) %in% rois2follow$model,"bold","plain"),
                       size=if_else(levels(d_csa$model) %in% rois2follow$model,11,10))
  ) +
  ## content - betas
  geom_point(size=3.5,alpha=0.8) +
  geom_errorbar(aes(xmin=Estimate-1.96*StdError, xmax=Estimate+1.96*StdError,width=.5), position=position_dodge(.0)) +
  ## more general styling
  scale_color_manual(values=brewer.pal(12, "Paired")[c(4,3)],
                     labels=c(
                       bquote("Model1a: Reading  ~ covariates + " ~ CSA[RH]),
                       bquote("Model2a: Reading  ~ covariates + " ~ totalCSA[RH] + CSA[RH]) )
  ) +
  scale_shape_manual(values=c(19,17),
                     labels=c(
                       bquote("Model1a: Reading  ~ covariates + " ~ CSA[RH]),
                       bquote("Model2a: Reading  ~ covariates + " ~ totalCSA[RH] + CSA[RH]) )
  ) +
  facet_grid(. ~ adjCognitive) +
  labs(title="",y="",color="Model",shape="Model") +
  theme(legend.position="bottom",legend.direction = "vertical") +
  # same x range for csa and ct
  coord_cartesian(xlim=c(min(d$Estimate)-0.02,max(d$Estimate)+0.02) %>% round(digits=2)) +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  theme(axis.title.y=element_blank(),legend.position="bottom",legend.direction = "vertical") +
  NULL

ct_betas <- ggplot(data=d_ct,
                   aes(x=Estimate,y=region,color=adjGlobal,shape=adjGlobal)) + 
  geom_vline(xintercept = 0,color="black",alpha=0.7,linetype="solid") +
  theme(axis.text.y=
          element_text(face=if_else(levels(d_ct$model) %in% rois2follow$model,"bold","plain"),
                       size=if_else(levels(d_ct$model) %in% rois2follow$model,11,10))) +
  ## content - betas
  geom_point(size=4,alpha=0.8) +
  geom_errorbar(aes(xmin=Estimate-1.96*StdError, xmax=Estimate+1.96*StdError,width=.5), position=position_dodge(.0)) +
  ## more general styling
  scale_color_manual(values=brewer.pal(12, "Paired")[c(8,7)],
                     labels=c(
                       bquote("Model1a: Reading  ~ covariates + " ~ CT[RH]),
                       bquote("Model2a: Reading  ~ covariates + " ~ meanCT[RH] + CT[RH]) )
                       ) +
  scale_shape_manual(values=c(19,17),
                     labels=c(
                       bquote("Model1a: Reading  ~ covariates + " ~ CT[RH]),
                       bquote("Model2a: Reading  ~ covariates + " ~ meanCT[RH] + CT[RH]) )
  ) +
  facet_grid(. ~ adjCognitive) +
  labs(title="",y="",color="Model",shape="Model") +
  # same x range for csa and ct
  coord_cartesian(xlim=c(min(d$Estimate)-0.02,max(d$Estimate)+0.02) %>% round(digits=2)) +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  theme(axis.title.y=element_blank(),legend.position="bottom",legend.direction = "vertical") +
  NULL
  

all_betas<-plot_grid(
  csa_betas,
  ct_betas,
  rel_heights = c(1.5,3),
  ncol=1, labels=c("A","B"))

# combined plots
out_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/","/output/")
if(!dir.exists(out_dir)){dir.create(out_dir)}

ggsave(all_betas,file=paste(out_dir,pheno,"_",atlas,"_rois2follow_sensitivity",".pdf",sep=""),height=7,width=13)


