# clean workspace
rm(list=ls())
# libraries and custom functions
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(RColorBrewer)
library(dplyr)
library(ggsignif)
#---------------------------------------------------
options(stringsAsFactors = FALSE)
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")
ldsc_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/ldsc/",sep="")
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
# read ggseg images and color table, for consistency
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/",sep="")
load(paste(pheno_dir,"ggseg_brain_ROIs.Rdata",sep=""))
braincolor_codes <- braincolor_codes %>% 
  mutate(region_name=Region  %>% gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% toupper() )



brains<-plot_grid(
  gg_area_all,
  gg_legend_thick_all,
  # plot_grid(titleA,gg_area_all,ncol=1,rel_heights = c(0.1,1),align="v"),
  # plot_grid(titleT,gg_thick_all,ncol=1,rel_heights = c(0.1,1)),
  nrow=1,rel_widths = c(1,1)
)
brains2<-plot_grid(
  gg_area_all,
  gg_legend_thick_all,
  # plot_grid(titleA,gg_area_all,ncol=1,rel_heights = c(0.1,1),align="v"),
  # plot_grid(titleT,gg_thick_all,ncol=1,rel_heights = c(0.1,1)),
  ncol=1 #,rel_widths = c(1,1)
)


#----------------------------------------------------------------------

#----------------------------------------------------------------------
# read data
h2<-read.csv(paste(ldsc_dir,"ldsc_h2_UKB_BIG40.csv",sep="")) %>% 
  mutate(region=factor(region,levels=unique(braincolor_codes$region_name))) %>% 
  arrange(region) %>% filter(hemisphere=="L")
rg<-read.csv(paste(ldsc_dir,"ldsc_rgLR_UKB_BIG40.csv",sep=""))  %>% 
  mutate(region=factor(region,levels=unique(braincolor_codes$region_name))) %>% arrange(region)
rg_cogn<-read.csv(paste(ldsc_dir,"ldsc_rg_cognitive_UKB_BIG40.csv",sep=""))  %>% 
  mutate(region=factor(region,levels=unique(braincolor_codes$region_name))) %>% arrange(region) %>%
  filter(!is.na(region)) %>% filter(hemisphere=="L")
# add labels for significance
h2$sig<-""
h2$sig[h2$p<0.05]<-"*"
h2$sig[h2$p<0.05/NROW(h2)]<-"***"

rg$sig<-""
rg$sig[rg$p2<0.05]<-"*"
rg$sig[rg$p2<0.05/NROW(rg)]<-"***"

rg_cogn$sig<-""
rg_cogn$sig[rg_cogn$p<0.05]<-"*"
rg_cogn$sig[rg_cogn$p<0.05/NROW(rg_cogn)]<-"***"
rg_cogn$cognitive<-factor(rg_cogn$cognitive, 
                          levels=c("Educational Attainement","Cognitive Performance","Developmental Dyslexia"))

rg_cogn$measure<-rg_cogn$measure %>% tolower() %>% sapply(.,simpleCap)
rg_cogn<-rg_cogn %>% arrange(cognitive,measure,region)
#----------------------------------------------------------------------
## def colors
mytheme<-  theme(axis.title.x=element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1,face="bold"),
                 axis.text=element_text(size=12),
                 legend.position="bottom")
rois<-subset(h2,!is.na(color))
# rois
colors_a <- rois  %>% filter(hemisphere=="L" & measure=="AREA") %>% pull(color)
colors_t <- rois  %>% filter(hemisphere=="L" & measure=="THICKNESS") %>% pull(color)
colors <- rois %>% pull(color) %>% unique()

# # add color tag to the tables - for consistency across plots
# for (d in c("h2","rg","rg_cogn","rg_read")){
#   dt<-get(d)
#   dt<-merge(dt,braincolor_codes)
#   # specify order of levels for consistency of colors, and order in plots
#   dt$region<-factor(dt$region,
#                     levels=toupper(levels(braincolor_codes$Region)))
#   dt <- dt %>% arrange(region) %>% mutate(region_name=factor(toupper(region),levels=unique(region)))
#   assign(d,dt)
# }
#----------------------------------------------------------------------
# edit names
h2$measure<-h2$measure %>% tolower() %>% sapply(.,simpleCap)
rg$measure<-rg$measure %>% tolower() %>% sapply(.,simpleCap)

#----------------------------------------------------------------------
# plots
rg_cogn_plot2 <- rg_cogn %>%
  ggplot(., aes(x=region,y=rg,width=0.7,color=region,fill=region)) +
  # geom_hline(yintercept =0,color="black",alpha=0.7,linetype="solid") +
  
  # annotate area and thickness, and add line in between
  annotate("text", x=levels(rg_cogn$region)[3],y=0.5,label="CSA",size=4) +
  annotate("text", x=levels(rg_cogn$region)[4],y=0.5,label="CT",size=4) +
  geom_vline(xintercept = 7-0.5,color="grey85",linetype="solid",size=1) +
  geom_hline(yintercept = 0,color="black",linetype="solid",size=1,alpha=0.6) +
  
  # plot data
  geom_errorbar(aes(ymin=rg.CI95upper, ymax=rg.CI95lower,width=.0), position=position_dodge(.7),stat="identity",size=1.3) +
  geom_point(position=position_dodge(.7), stat="identity",size=4,shape=24,color="lightgrey") +
  geom_text(aes(x=region,y=rg.CI95upper+0.1,label=sig),position=position_dodge(width=0.7),size=10, stat="identity") +
  
  # legend and axes
  # ylab(bquote('Estimate ('*rho*')')) +
  ylab(bquote(r[g])) +
  scale_color_manual(values =colors,guide="none") +
  scale_fill_manual(values =colors,guide="none") +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_discrete(limits=rev(levels(rg_cogn$region))) +
  scale_y_continuous(breaks=seq(-0.6,0.6,by=0.30)) +
  # 
  coord_flip(ylim=c(-0.6,0.6)) +
  facet_grid( .~cognitive, switch="y") +
  # theme adjustements
  theme_minimal_vgrid(16,line_size=0.3) +
  panel_border() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        # facets
        strip.background = element_blank(),
        # remove legend
        legend.position="none") +
  NULL


rg_cogn_plot2
brain_rgL_cogn_plot22<-plot_grid(brains,rg_cogn_plot2,nrow=2,rel_heights = c(1.3,2))
ggsave(brain_rgL_cogn_plot22,file=paste0(out_dir,"brainL_rg_leftROIs_cogn_ukb_big40_slides.png"),
       width=10,height=8)
