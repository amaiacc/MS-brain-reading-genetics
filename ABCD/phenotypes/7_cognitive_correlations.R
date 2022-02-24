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
library(dplyr); library(tidyr)
library(data.table)
library(psych)
#
library(scales)
library(ggplot2)
# if (!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("caijun/ggcorrplot2")
# library(ggcorrplot2)
library(ggcorrplot)
library(corrplot)
library(cowplot)

#
corrtheme <- theme(
  axis.title.x = element_text(angle = 0, color = 'grey20',face="bold", size=12),
  axis.title.y = element_text(angle = 90, color = 'grey20',face="bold", size=12),
  legend.title=element_text(face="bold", size=12),
  plot.title = element_text(angle = 0, face="bold", size=14),
  axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 10),
  axis.text.y = element_text(angle=0, hjust = 1,vjust=1,size = 10)
) 
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
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",sep="")
input_dir=paste(pheno_dir,"baseline","_",atlas,"/",sep="") %>% gsub(" ","",.)
working_dir=paste(pheno_dir,atlas,sep="")
if(!dir.exists(working_dir)){dir.create(working_dir)}
setwd(working_dir)
#---------------------------------------------------

#---------------------------------------------------
# read data
f=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.csv",sep="")
d<-read.csv(f)
d$rel_family_id <- as.factor(d$rel_family_id)
colnames(d) <- colnames(d) %>% gsub("smri_|cort.destrieux_","",.) %>% gsub("\\.|_","_",.)
# residuals
f2=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis_ROIs2follow_residuals.csv",sep="")
d_res<-read.csv(f2)
colnames(d_res) <- colnames(d_res) %>% gsub("smri_|cort.destrieux_","",.) %>% gsub("\\.|_","_",.)

rm(f,f2)
#---------------------------------------------------
# define cognitive variables
vars2follow <- colnames(d)[grep("nihtbx|wiscv",colnames(d))]
#---------------------------------------------------

#---------------------------------------------------
# correlations across cognitive variables + visualize
#---------------------------------------------------


#---------------------------------------------------
# RAW measures
#---------------------------------------------------
d1<-d[,vars2follow]
colnames(d1) <- colnames(d1) 
corROIS<-corr.test(d1,minlength =100,adjust="bonferroni")

## to long format
corROIs_long<-do.call("cbind",
  lapply(c("r","se","p"), function(m){
    x <- corROIS[m] %>% as.data.frame()
    x2 <- x %>%
      mutate(measure1=rownames(x)) %>%
      pivot_longer(cols=colnames(x),names_to="measure2",values_to=m)
    return(x2)
    } )
  ) %>% subset(select=which(!duplicated(names(.)))) %>% 
  mutate(
    measure2=gsub("^r.","",measure2),
    n=corROIS$n, adj=corROIS$adjust,type="raw measures") %>%
  # remove rows of phenotypes r with itself
  filter(measure1!=measure2) %>%
  distinct()

tmp1<- corROIs_long %>% rename(m1=measure1, m2=measure2)
tmp2<- corROIs_long %>% rename(m1=measure2, m2=measure1)
tmp<-rbind(tmp1,tmp2)
tmp[which(duplicated(tmp)),] %>% View()
## plot
c2<-corROIS$r
p2<-corROIS$p

# change names
names<-colnames(c2) %>% gsub(".lh|.rh|uncorrected|tss","",.) %>% gsub("\\.|_"," ",.) %>%
  gsub("nihtbx ","",.) %>% gsub(" $","",.) %>%
  sapply(simpleCap)
names<-names %>% gsub("Picvocab","Vocabulary",.) %>%
  gsub("Fluidcomp","FluidComponent",.) %>% 
  gsub("Fluidcomp","FluidComponent",.) %>%
  gsub("Pea Wiscv","WISCV",.)
colnames(c2)<-rownames(c2)<-names
colnames(p2)<-rownames(p2)<-names
colnames(d1)<-names
rm(names)

names<-c("Reading","Vocabulary","FluidComponent","WISCV")
c2<-c2[names,names]
p2<-p2[names,names]
corrplot_unadj <- ggcorrplot::ggcorrplot(c2,
           p.mat=p2,
           method="square",
           insig="blank",
           hc.order=FALSE,
           lab=TRUE,
           type="lower",
           show.diag=FALSE,
           digits=2
           # colors = c("#6D9EC1", "white", "#E46726"),
           ) + 
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"),
                       midpoint = 0,
                       limits=c(-1,1)) +
  corrtheme +
  labs(y = '', 
       x = '',
       fill="Corr. Coef.") + 
  NULL

library(hexbin)
ggpairs_hex <- function(df, hexbins = 10) {
  # REF: https://stackoverflow.com/questions/20872133/using-stat-binhex-with-ggpairs
  p <- ggpairs(df, lower="blank")
  seq <- 1:ncol(df)
  for (x in seq)
    for (y in seq) 
      if (y>x) 
        p <- putPlot(p, ggplot(df, aes_string(x=names(df)[x],y=names(df)[y])) + stat_binhex(bins = hexbins), y,x)
  
  return(p)
}

corrpairs <- d1[,names] %>% ggpairs_hex(., hexbins=20) + theme_bw()
#---------------------------------------------------

out_dir=paste0(working_dir,"/summary/figures/")

corrplot_unadj %>%
  ggsave(file=paste0(out_dir,"corrplot_cognitive_all_unadj.pdf"))

corrpairs %>% ggsave(file=paste0(out_dir,"corrpairs_cognitive_all_unadj.png"),height=4.5,width=4.5)

