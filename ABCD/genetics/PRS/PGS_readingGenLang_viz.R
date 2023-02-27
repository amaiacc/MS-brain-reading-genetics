# ** to add data etc, continuing after PGS_lmm.R **
rm(list=ls())
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects/resources/datasets/ABCD/scripts/genetics/PRS/base_readingGenLang.config"


# plotting
library(ggplot2)
library(cowplot);  theme_set(theme_cowplot())
library(ggridges)
library(psych)
library(RColorBrewer)
# data arranging
library(tidyr)
library(dplyr)

# matSpD:
matSpD_file="F:/projects/resources/matSpD.R"
# set options
options(stringsAsFactors=FALSE)

#---------------------------------------------------
# QC PRS
#---------------------------------------------------
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/resources/datasets/"} else {dir="/export/home/acarrion/acarrion/projects/resources/datasets/"}
# source config file to define parameters for this run
source(config_file)
td=paste(project,"_",batch,sep="")

# define paths to data
prs_dir=paste(dir,project,"/data/working_data/genotyping/prsice/",base_project,"/",sep="")
out_dir=paste(prs_dir,"/output/",sep="")
if(!dir.exists(out_dir)){dir.create(out_dir)}
# set working dir
setwd(prs_dir)

#---------------------------------------------------
# read data
#---------------------------------------------------
# pgs + genetic PCs (for all and EUR only)
pgs_summary<- read.csv(paste(out_dir,"PGS_",base_project,"_",s,"_summary_stats.csv",sep=""))
pgs_summary$threshold<- pgs_summary$threshold %>% as.character %>% as.numeric %>% as.factor()
pgs_summary$pheno<-pgs_summary$target_pheno %>% gsub("_","\n",.) %>% gsub("\\.lh","\nlh",.) %>% gsub("\\.rh","\nrh",.)
pgs_summary$measure<-NA
pgs_summary$measure[grep("area",pgs_summary$phenotype)]<-"area"
pgs_summary$measure[grep("thick",pgs_summary$phenotype)]<-"thickness"


pgs_summary$adj[grep("+ smri_vol_cort.destrieux_total",pgs_summary$base_model)]<-"tBV"
pgs_summary$adj[grep("+ smri_area_cort.destrieux_total",pgs_summary$base_model)]<-"tCSA"
pgs_summary$adj[grep("+ smri_thick_cort.destrieux_mean",pgs_summary$base_model)]<-"mCT"

# define brain and cognitive/behavioural measures
vars2test_c<-pgs_summary$phenotype[grep("nihtbx",pgs_summary$phenotype)] %>% unique()

vars2test_b<-pgs_summary$phenotype[grep("smri",pgs_summary$phenotype)] %>% unique()
vars2test_b_global<-vars2test_b[grep("total|mean",vars2test_b)] %>% unique()
vars2test_b_rois<-vars2test_b[grep("total|mean",vars2test_b,invert=TRUE)] %>% unique()

  
#------------------------------------
# visualize
#------------------------------------
## READING on BEHAVIOUR (cognition)
  
readRT_prs_plot_c<-ggplot(data=subset(pgs_summary,base_pheno=="WR.RT"&subset==s) %>%  filter(phenotype %in% vars2test_c),
                      aes(x=threshold,y=R2*100,fill=-log10(P))) +
  geom_col(width=1,position="dodge") + 
  # facet_grid(.~pheno,scale="free") +
  scale_fill_distiller(palette = 3, direction=1, name=expression(-log[10]*"(P-value)"),  
                       n.breaks=10,
                       guide = guide_legend(direction = "vertical", title.position = "right", 
                                            title.theme = element_text(angle = 90),title.vjust=0.5) ) +
  # scale_fill_distiller(palette = 9, direction=1, name=expression(-log[10]*"(P-value)"),  
  #                      guide = guide_legend(direction = "vertical", title.position = "right", title.theme = element_text(angle = 90),title.vjust=0.5) ) +
  labs(y=expression("% variance " (Delta~R^2)),x="GWAS P-value threshold",
       title="PGS Reading RT") +
  # mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL


# readACC_prs_plot_c<-ggplot(data=subset(pgs_summary,base_pheno=="WR.Z"&subset==s) %>%  filter(phenotype %in% vars2test_c),
#                       aes(x=threshold,y=R2*100,fill=-log10(P))) +
#   geom_col(width=1,position="dodge") + 
#   facet_grid(.~target_pheno,scale="free") +
#   # scale_fill_distiller(palette = 9, direction=1, name=expression(-log[10]*"(P-value)"),  
#   #                      guide = guide_legend(direction = "vertical", title.position = "right", title.theme = element_text(angle = 90),title.vjust=0.5) ) +
#   # labs(y=expression("% variance " (Delta~R^2)),x="GWAS P-value threshold",
#        title="Reading accuracy (GenLang GWASMA)") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   NULL

# prs_all_plot_c<-plot_grid(readACC_prs_plot_c,readRT_prs_plot_c,nrow=2)
prs_all_plot_c<-readRT_prs_plot_c
#---------------------------------------------------
# save summary plot and table
#---------------------------------------------------
ggsave(prs_all_plot_c,file=paste0(out_dir,"PRS_effects_readingGenlang2ABCD_",s,"_beh.pdf"),width=5,height=3)

