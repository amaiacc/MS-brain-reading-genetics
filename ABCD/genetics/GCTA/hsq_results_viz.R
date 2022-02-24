rm(list=ls())
library(tidyr)
library(dplyr)
library(tidyverse)
# library(ggpubr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(reshape2)
library(grid)
library(gridExtra)
library("ggseg")
library("ggsegExtra")
# ggseg::install_atlases()
library(ggsegDesterieux)
library(ggsegDefaultExtra)
#----------------------------------------------------------------------
#------------------------------------------
# plots, some general parameters
#------------------------------------------
options(stringsAsFactors = FALSE)
# define pattern for this run
root="ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned_rm05_adj"

# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
scripts_dir=paste(dir,"resources/datasets/ABCD/scripts/genetics/GCTA/",sep="")
working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/GCTA/",sep="")
out_dir=paste(working_dir,"/output/",sep="")
#
dir.create(file.path(out_dir), showWarnings = FALSE)
# out_dir="K://written/tex/multilateral/presentations/figures/"
setwd(out_dir)
source(paste(scripts_dir,"helper_functions.R",sep=""))
# read ggseg images and color table, for consistency
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/",sep="")
load(paste(pheno_dir,"ggseg_brain_ROIs.Rdata",sep=""))
braincolor_codes<- braincolor_codes %>% mutate(region_name=Region  %>% gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% toupper() )

#----------------------------------------------------------------------


#----------------------------------------------------------------------
h2<-read.csv(paste(out_dir,"gcta_estimates_","h2","_",root,".csv",sep=""))
h2$estimate<-h2$h2
rg<-read.csv(paste(out_dir,"gcta_estimates_","rg","_",root,".csv",sep=""))
rg$estimate<-rg$rg
#
# add labels for significance
h2$sig<-""
h2$sig[h2$pval<0.05]<-"*"
h2$sig[h2$pval<0.05/NROW(h2)]<-"***"
rg$sig<-""
rg$sig[rg$pval<0.05]<-"*"
rg$sig[rg$pval<0.05/NROW(rg)]<-"***"
#----------------------------------------------------------------------
summary<-merge(h2,rg,all=TRUE) %>% select(-h2,-rg)
summary<-subset(summary,run=="residuals")

# define categories for the phenotypes
summary$category<-NA
summary$category[grep("smri*.*total",summary$p1)]<-"Brain global"
summary$category[grep("smri",summary$p1)[grep("total",summary$p1[grep("smri",summary$p1)],invert=TRUE)]]<-"Brain ROI"
summary$category[grep("smri",summary$p1,invert=TRUE)]<-"Cognitive"

# clean phenotype names
summary$pheno1<-summary$p1 %>% gsub("nihtbx_|pea_|_tss|_uncorrected|smri_|_cort.|destrieux|_adjGlobal","",.)
summary$pheno2<-summary$p2 %>% gsub("nihtbx_pea_|_tss|_uncorrected|smri_|_cort.|destrieux|_adjGlobal","",.)

# clean region name
summary$region<-gsub(".rh|.lh","",summary$pheno1)
summary$region_name <- summary$region %>% gsub("smri_thick_cort.|smri_area_cort.|destrieux_","",.) %>% 
  gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% gsub("AREA_|area_|THICK_|thick_","",.) %>% toupper()
summary$measure[grep("area",summary$pheno1)]<-"AREA"
summary$measure[grep("thick",summary$pheno1)]<-"THICKNESS"
# 
summary<-merge(summary,braincolor_codes%>% select(-hemisphere),by=c("region_name","measure"),all.x=TRUE) 
#----------------------------------------------------------------------
# PLOTS
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# plot theme
mytheme<-  theme(axis.title.x=element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1,face="bold"),
                 axis.text=element_text(size=12),
                 legend.position="bottom")
## ROIs and colors
rois <- subset(summary,(category=="Brain ROI"|category=="Brain global")&region_name!="VOLDESIKAN_TOTAL") %>% 
  mutate(region_name=factor(region_name,levels=unique(braincolor_codes$region_name))) %>% 
  arrange(region_name)

colors_a <- rois %>% filter(stat=="h2") %>% filter(hemisphere=="L" & measure=="AREA" & adj=="-") %>% pull(color)
colors_t <- rois %>% filter(stat=="h2") %>% filter(hemisphere=="L" & measure=="THICKNESS" & adj=="-") %>% pull(color)
#----------------------------------------------------------------------
h2_plot_roisLarea <- ggplot(data=rois %>% filter(stat=="h2") %>% 
                              filter(hemisphere=="L"&measure=="AREA"),
                            aes(x=region_name,y=estimate,width=0.7,
                                color=region_name,fill=region_name)) + 
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity",alpha=0.5,size=1) +
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region_name,y=estimate+1.96*se+0.1,label=sig),position=position_dodge(width=0.7),size=10, stat="identity") +
  ylab(bquote('Estimate ('*h^2*')')) + xlab("") + 
  facet_grid(.~ measure,scales = "free_x",space="free_x") +
  scale_color_manual(values = colors_a,guide="none") +
  scale_fill_manual(values = colors_a,guide="none") +
  facet_grid( adj  ~ .) +
  # coord_flip() +
  ylim(c(-0.1,1)) +
  mytheme +
  NULL

h2_plot_roisLthick <- ggplot(data=rois %>% filter(stat=="h2") %>% filter(hemisphere=="L"&measure=="THICKNESS"),
                           aes(x=region_name,y=estimate,width=0.7,
                               color=region_name,fill=region_name)) + 
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity",alpha=0.5,size=1) +
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region_name,y=estimate+1.96*se+0.1,label=sig),position=position_dodge(width=0.7),size=10, stat="identity") +
  ylab(bquote('Estimate ('*h^2*')')) + xlab("") + 
  facet_grid(.~ measure,scales = "free_x",space="free_x") +
  scale_color_manual(values =colors_t,guide="none") +
  scale_fill_manual(values =colors_t,guide="none") +
  facet_grid( adj  ~ .) +
  # coord_flip() +
  ylim(c(-0.1,1)) +
  mytheme +
  NULL


row1<-plot_grid(gg_area_all,NULL,gg_thick_all,nrow=1,rel_widths = c(2,0.5,2))
row2<-plot_grid(NULL,h2_plot_roisLarea,NULL,h2_plot_roisLthick,
                rel_widths = c(0.1,1.2,0.5,2),
                align="h",nrow=1)
gg_h2_all<-plot_grid(row1,NULL,row2,ncol=1,rel_heights = c(1.2,0.05,2))


# heritabilities
colors<-rois %>% filter(stat=="h2" & !is.na(hemisphere)) %>% pull(color) %>% unique()
h2_plot_rois<-ggplot(data=rois %>% filter(stat=="h2" & !is.na(hemisphere)),
                             aes(x=region_name,y=estimate,width=0.7,
                                 color=region_name,shape=hemisphere)) + 
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_bar(aes(y=estimate,alpha=hemisphere,fill=region_name),position=position_dodge(), stat="identity",size=1) +
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region_name,y=estimate+1.96*se+0.1,label=sig),position=position_dodge(width=0.7),size=10, stat="identity") +
  ylab(bquote('Estimate ('*h^2*')')) + xlab("") + 
  facet_grid(.~ measure,scales = "free_x",space="free_x") +
  scale_color_manual(values =colors,guide="none") +
  scale_fill_manual(values =colors,guide="none") +
  scale_alpha_discrete(range=c(0.1,0.8)) +
  facet_grid( adj  ~ .) +
  # coord_flip() +
  ylim(c(-0.1,1)) +
  mytheme +
  NULL

h2_plot_roisL<-ggplot(data=rois %>% filter(stat=="h2") %>% filter(hemisphere=="L"),
                     aes(x=region_name,y=estimate,width=0.7,
                         color=region_name,fill=region_name)) + 
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity",alpha=0.5,size=1) +
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region_name,y=estimate+1.96*se+0.1,label=sig),position=position_dodge(width=0.7),size=10, stat="identity") +
  ylab(bquote('Estimate ('*h^2*')')) + xlab("") + 
  facet_grid(.~ measure,scales = "free_x",space="free_x") +
  scale_color_manual(values =colors,guide="none") +
  scale_fill_manual(values =colors,guide="none") +
  facet_grid( adj  ~ .) +
  # coord_flip() +
  ylim(c(-0.1,1)) +
  mytheme +
  NULL


# genetic correlations between left and right
rg_plot_rois <- ggplot(data=rois %>% filter(stat=="rg"),
                  aes(x=region_name,y=estimate,width=0.7,color=region_name)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="solid") +
  geom_hline(yintercept = 1,color="black",alpha=0.7,linetype="dashed") +
  geom_hline(yintercept = -1,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3) +
  geom_bar(aes(y=estimate),position=position_dodge(0.3), stat="identity",fill="white",size=1,width = 0.1) +
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*rho*')')) +
  facet_grid(.~ measure,scales = "free_x",space="free_x") +
  scale_color_manual(values =colors,guide="none") +
  # coord_flip() +
  facet_grid( adj  ~ .) +
  coord_cartesian(ylim=c(-1.1,1.1)) +
  mytheme +
  NULL

## cognitive
# heritabilities
h2_plot_cogn <- ggplot(data=subset(summary,stat=="h2"&category=="Cognitive") %>%
                         mutate(region_name=factor(region_name,
                                                      levels=c("READING","PICVOCAB","CRYST","FLUIDCOMP","WISCV"))),
                     aes(x=region_name,y=estimate,width=0.7,
                         color=region_name,fill=region_name)) + 
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity",size=1,alpha=0.8) +
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region_name,y=estimate+1.96*se+0.1,label=sig),position=position_dodge(width=0.7),size=10, stat="identity") +
  ylab(bquote('Estimate ('*h^2*')')) + xlab("") + 
  scale_color_brewer(palette="Set2",guide="none") +
  scale_fill_brewer(palette="Set2",guide="none") +
  # coord_flip() +
  ylim(c(-0.1,1)) +
  mytheme +
  NULL

#----------------------------------------------------------------------
# save plots
## h2 + brains
row12<-plot_grid(NULL,gg_area_all,NULL,gg_thick_all,NULL,nrow=1,rel_widths = c(0.3,2,0.3,2,0.3))
ggsave(plot_grid(row12, NULL,h2_plot_rois,ncol=1,rel_heights = c(1.2,0.1,2)),file="h2_ROIs_gg_abcd.png")
ggsave(plot_grid(row12, NULL,h2_plot_roisL,ncol=1,rel_heights = c(1.2,0.1,2)),file="h2_leftROIs_gg_abcd.png")
ggsave(gg_h2_all,file="h2_leftROIs_gg_abcd_v2.png")
# h2 and rg
ggsave(h2_plot_roisL,file="h2_leftROIs_abcd.png")
ggsave(h2_plot_rois,file="h2_ROIs_abcd.png")
ggsave(rg_plot_rois,file="rgLR_ROIs_abcd.png")
# cognitive
ggsave(h2_plot_cogn,file="h2_cognitive_abcd.png")
#----------------------------------------------------------------------




# # build corrplots
# library(corrplot)
# 
# ## across cognitive variables
# # complete matrix
# tmp1<- rg %>% filter(str_detect(p1,"nihtbx") & str_detect(p2,"nihtbx") ) %>% select(p1,p2,rg,se,pval)
# tmp2<- tmp1 %>% select(p2,p1,rg,se,pval)
# colnames(tmp2)<-colnames(tmp1)
# tmp<-rbind(tmp1,tmp2)
# rm(tmp1,tmp2)
# # to wide
# t_rg<-tmp %>% select(p1,p2,rg) %>% spread(p2,rg) %>% select(-1) %>% as.matrix()
# t_p<-tmp %>%  select(p1,p2,pval) %>% spread(p2,pval) %>% select(-1) %>% as.matrix()
# rownames(t_rg)<-rownames(t_p)<-colnames(t_rg)
# 
# corrplot(corr=t_rg, p.mat=t_p,
#          method="circle",
#          order="hclust",
#          # style
#          col=brewer.pal(n=8, name="PuOr"),
#          insig="blank",
#          is.corr=FALSE,
#          addCoef.col=TRUE,
#          addCoefasPercent=TRUE,
#          type="lower",
#          diag=FALSE,
#          tl.col="black",
#          tl.srt=45,
#          sig.level=0.05,
#          mar=c(0,0,1,0),
#          title="Cognitive")
# 
# ## across brain variables
# 
# # complete matrix
# tmp1<- rg %>% filter(str_detect(p1,"smri") & str_detect(p2,"smri") ) %>% select(p1,p2,rg,se,pval)
# tmp2<- tmp1 %>% select(p2,p1,rg,se,pval)
# colnames(tmp2)<-colnames(tmp1)
# tmp<-rbind(tmp1,tmp2)
# rm(tmp1,tmp2)
# # to wide
# w<-which(duplicated(tmp %>% select(p1,p2,rg) ))
# tmp<-tmp[-w,]
# t_rg<-tmp[-w,] %>% select(p1,p2,rg) %>% spread(p2,rg) %>% select(-1) %>% as.matrix()
# t_p<-tmp[-w,] %>%  select(p1,p2,pval) %>% spread(p2,pval) %>% select(-1) %>% as.matrix()
# rownames(t_rg)<-rownames(t_p)<-colnames(t_rg)
# 
# ## now make lower triangle ->LH and upper triangel ->RH
# 
# 
# corrplot(corr=t_rg, p.mat=t_p,
#          method="circle",
#          order="hclust",
#          # style
#          col=brewer.pal(n=8, name="PuOr"),
#          insig="label_sig",
#          is.corr=TRUE,
#          # addCoef.col=TRUE,
#          # addCoefasPercent=TRUE,
#          type="lower",
#          diag=FALSE,
#          tl.col="black",
#          tl.srt=45,
#          sig.level=0.05,
#          mar=c(0,0,1,0),
#          title="Cognitive")
# 
