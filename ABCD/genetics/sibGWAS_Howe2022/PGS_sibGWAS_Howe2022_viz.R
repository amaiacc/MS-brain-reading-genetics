# ** to add data etc, continuing after PGS_lmm.R **
rm(list=ls())
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects/resources/datasets/ABCD/scripts/genetics/sibGWAS_Howe2022/base_sibGWAS_Howe2022.config"

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
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
# QC PGS
#---------------------------------------------------
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/resources/datasets/"} else {dir="/export/home/acarrion/acarrion/projects/resources/datasets/"}
# source config file to define parameters for this run
source(config_file)
td=paste(project,"_",batch,sep="")

source(paste0(dir,"/../../general_scripts/helper_functions.R"))

# define paths to data
PGS_dir=paste(dir,project,"/data/working_data/genotyping/prsice/",base_project,"/",sep="")
out_dir=paste(PGS_dir,"/output/",sep="")
if(!dir.exists(out_dir)){dir.create(out_dir)}
# set working dir
setwd(PGS_dir)

# read ggseg images and color table, for consistency
pheno_dir=paste(dir,"/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/",sep="")
load(paste(pheno_dir,"ggseg_brain_ROIs.Rdata",sep=""))
# rois colors
colors_a <- braincolor_codes  %>% filter(hemisphere=="L" & measure=="AREA") %>% pull(color)
colors_t <- braincolor_codes  %>% filter(hemisphere=="L" & measure=="THICKNESS") %>% pull(color)
colors<- braincolor_codes %>% pull(color) %>% unique()  

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

# list of base phenotypes
howe_dir=paste0(dir,"/GWAS_sumstats/downloaded_data/ieu-openGWAS/")
howe2022_traits<-read.table(paste0(howe_dir,"Howe2022_GWASlist2download.txt"),header=TRUE) %>%
  mutate(study=paste0(strsplit(author," ") %>% sapply("[[",1),year),
         measure.Howe2022=paste(trait,note),
         id=gsub("-","\\.",id)
  ) %>%
  dplyr::select(study,id,trait,note,measure.Howe2022,sample_size)
#
pgs_summary<- merge(pgs_summary,howe2022_traits,by.x="base_pheno",by.y="id",all.x = TRUE)


#------------------------------------
# visualize
#------------------------------------
## COGN on BEHAVIOUR
max_p<- (-log10(subset(pgs_summary,subset==s)$P)) %>% max() %>% round(digits=0) +1
base_p="Years of schooling"
ys_PGS_plot_c<-ggplot(data=subset(pgs_summary,trait==base_p&subset==s) %>%  
                        filter(phenotype %in% vars2test_c),
                      aes(x=threshold,y=R2*100,fill=-log10(P))) +
  geom_col(width=1,position="dodge") + 
  facet_grid(.~note,scale="free") +
  scale_fill_distiller(palette = 3, direction=1, name=expression(-log[10]*"(P-value)"),
                       limit=c(0,max_p),n.breaks=10,
                       guide = guide_legend(direction = "vertical", title.position = "right", 
                                            title.theme = element_text(angle = 90),title.vjust=0.5) ) +
  labs(y=expression("% variance " (Delta~R^2)),x="GWAS P-value threshold",
       title=paste0("PGS ",base_p)) +
  # mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim=c(0,1)) +
  NULL

base_p="Cognitive function"
cf_PGS_plot_c<- ggplot(data=subset(pgs_summary,trait==base_p&subset==s) %>%  
                         filter(phenotype %in% vars2test_c),
                       aes(x=threshold,y=R2*100,fill=-log10(P))) +
  geom_col(width=1,position="dodge") + 
  facet_grid(.~note,scale="free") +
  scale_fill_distiller(palette = 3, direction=1, name=expression(-log[10]*"(P-value)"),  
                       limit=c(0,max_p),n.breaks=10,
                       guide = guide_legend(direction = "vertical", title.position = "right", 
                                            title.theme = element_text(angle = 90),title.vjust=0.5) ) +
  labs(y=expression("% variance " (Delta~R^2)),x="GWAS P-value threshold",
       title=paste0("PGS ",base_p)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim=c(0,1)) +
  NULL


thrlist=levels(pgs_summary$threshold)


## define best base dataset and p-value threshold for predicting reading
reading_best<-subset(pgs_summary,subset==s) %>%  filter(phenotype %in% vars2test_c)
reading_best<-reading_best %>% filter(R2==max(reading_best$R2))

best_base<-reading_best$base_pheno
best_thr<-reading_best$threshold
best_note<-reading_best$note

ys_PGS_plot_c2 <- ys_PGS_plot_c + geom_rect(aes(
  xmin=which(thrlist==best_thr)-0.5,
  xmax=which(thrlist==best_thr)+0.5,
  ymin =0, ymax = max(reading_best$R2)*100),
  color="black",size=1.5) +
  guides(fill=guide_legend(override.aes = list(color = "white"),
                           direction = "vertical", 
                           title.position = "right", 
                           title.theme = element_text(angle = 90),title.vjust=0.5)) +
  theme(axis.text.x=
          element_text(face=if_else(
            (levels(pgs_summary$threshold) %in% as.character(best_thr)) & (unique(pgs_summary$note) %in% best_note),
            "bold","plain"),
                       size=if_else(levels(pgs_summary$threshold) %in% as.character(best_thr),11,10))
  ) +
  NULL


PGS_legend<-get_legend(ys_PGS_plot_c)

PGS_all_plot_c<-plot_grid(ys_PGS_plot_c + theme(legend.position = "none"),
                          cf_PGS_plot_c + theme(legend.position = "none"),
                          NULL,
                          PGS_legend,
                          nrow=1,
                          labels=c("a","b",""),rel_widths = c(1,1,0.05,0.2))


#---------------------------------------------------
# save summary plot and table
#---------------------------------------------------
ggsave(PGS_all_plot_c,file=paste0(out_dir,"PGS_effects_",base_project,"_","2ABCD_",s,"_behav.pdf"),width=11,height=3)

# clean intermediate
rm(list=ls(pattern="plot"))

#---------------------------------------------------