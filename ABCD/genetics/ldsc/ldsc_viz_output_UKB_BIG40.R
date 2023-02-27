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
out_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/output/")
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
# read ggseg images and color table, for consistency
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/",sep="")
load(paste(pheno_dir,"ggseg_brain_ROIs.Rdata",sep=""))


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
rg_cogn1<-read.csv(paste(ldsc_dir,"ldsc_rg_cognitive_UKB_BIG40.csv",sep=""))  %>% 
  mutate(region=factor(region,levels=unique(braincolor_codes$region_name))) %>% arrange(region) %>%
  filter(!is.na(region)) %>% filter(hemisphere=="L")
rg_read<-read.csv(paste(ldsc_dir,"ldsc_rg_readingGenLang_UKB_BIG40.csv",sep=""))  %>% 
  mutate(region=factor(region,levels=unique(braincolor_codes$region_name))) %>% arrange(region) %>%
  filter(!is.na(region)) %>% filter(hemisphere=="L")

rg_sibs<-read.csv(paste(ldsc_dir,"ldsc_rg_sibGWASHowe2022_UKB_BIG40.csv",sep=""))  %>% 
  mutate(region=factor(region,levels=unique(braincolor_codes$region_name))) %>% arrange(region) %>%
  filter(!is.na(region)) %>% filter(hemisphere=="L")

rg_cogn<-merge(rg_cogn1,rg_read,all=TRUE)
rg_cogn<-merge(rg_cogn,rg_sibs,all=TRUE)

# create factor with cognitive trait
rg_cogn<-rg_cogn %>% 
  mutate(cognitive=if_else(is.na(cognitive),
                           gsub("_"," ",name) %>% gsub("WR","Word Reading",.) %>% gsub("RT","(RT)",.) %>% gsub("Z","(Accuracy)",.),
                           cognitive),
         cognitive=if_else(is.na(cognitive),
                           trait,
                           cognitive)
         ) %>%
  # mutate(cognitive=factor(cognitive,
  #                         levels=c("Educational Attainement",
  #                                  "Cognitive Performance",
  #                                  "Developmental Dyslexia",
  #                                  "Word Reading (RT)"))) %>%
  mutate(measure=tolower(measure) %>% sapply(.,simpleCap) )



# correct for multiple comparisons, FDR
(p.adjust(rg_cogn$p,method = "bonferroni")<0.05) %>% table()

# add labels for significance
h2$sig<-""
h2$sig[h2$p<0.05]<-"*"
h2$sig[h2$p<(0.05/NROW(h2))]<-"***"

rg$sig<-""
rg$sig[rg$p2<0.05]<-"*"
rg$sig[rg$p2<(0.05/NROW(rg))]<-"***"

rg_cogn$sig<-""
rg_cogn$sig[rg_cogn$p<0.05]<-"*"
rg_cogn$sig[rg_cogn$p<(0.05/NROW(rg_cogn))]<-"***"


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
rg_cogn_plot2 <- rg_cogn %>% filter(is.na(study)) %>%
  mutate(cognitive=factor(cognitive,
                          levels=c("Educational Attainement",
                                   "Cognitive Performance",
                                   "Developmental Dyslexia",
                                   "Word Reading (Accuracy)",
                                   "Word Reading (RT)"))) %>%
  filter(cognitive!="Word Reading (Accuracy)") %>%
  mutate(sig_val=if_else(rg<0,rg.CI95lower*1.2,rg.CI95upper*1.2)) %>%
  ggplot(., aes(x=region,y=rg,width=0.7,color=region,fill=region,shape=note,alpha=sig)) +
  # geom_hline(yintercept =0,color="black",alpha=0.7,linetype="solid") +
  
  # annotate area and thickness, and add line in between
  annotate("text", x=levels(rg_cogn$region)[3],y=0.55,label="CSA",size=4) +
  annotate("text", x=levels(rg_cogn$region)[4],y=0.55,label="CT",size=4) +
  geom_vline(xintercept = 7-0.5,color="grey85",linetype="solid",size=1) +
  geom_hline(yintercept = 0,color="black",linetype="solid",size=1,alpha=0.6) +
  
  # plot data
  geom_errorbar(aes(ymin=rg.CI95upper, ymax=rg.CI95lower,width=.0), position=position_dodge(.7),stat="identity",size=1.3) +
  geom_point(position=position_dodge(.7), stat="identity",size=4,shape=24,color="lightgrey") +
  geom_text(aes(x=region,y=sig_val,label=sig),position=position_dodge(width=0.7),size=10, stat="identity") +
  
  # legend and axes
  # ylab(bquote('Estimate ('*rho*')')) +
  ylab(bquote(r[g])) +
  scale_color_manual(values =colors,guide="none") +
  scale_fill_manual(values =colors,guide="none") +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_discrete(limits=rev(levels(rg_cogn$region))) +
  scale_y_continuous(breaks=seq(-0.6,0.6,by=0.30)) +
  # 
  coord_flip(ylim=c(-0.85,0.85)) +
  # coord_flip() +
  facet_grid( .~cognitive, switch="y",scales="free") +
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
ggsave(brain_rgL_cogn_plot22,file=paste0(out_dir,"brainL_rg_leftROIs_cognANDread_ukb_big40.png"),
       width=12,height=8)


rg_sibs_plot <- rg_cogn %>% filter(cognitive!="Height") %>% filter(!is.na(study)) %>%
  mutate(cognitive=factor(cognitive,
                          levels=c("Years of schooling",
                                   "Cognitive function",
                                   "Height"))) %>%
  mutate(sig_val=if_else(rg<0,rg.CI95lower*1.1,rg.CI95upper*1.1)) %>%
  ggplot(., aes(x=region,y=rg,width=0.7,color=region,fill=region,shape=note,alpha=sig)) +
  # geom_hline(yintercept =0,color="black",alpha=0.7,linetype="solid") +
  
  # annotate area and thickness, and add line in between
  annotate("text", x=levels(rg_cogn$region)[3],y=0.55,label="CSA",size=4) +
  annotate("text", x=levels(rg_cogn$region)[4],y=0.55,label="CT",size=4) +
  geom_vline(xintercept = 7-0.5,color="grey85",linetype="solid",size=1) +
  geom_hline(yintercept = 0,color="black",linetype="solid",size=1,alpha=0.6) +
  
  # plot data
  geom_errorbar(aes(ymin=rg.CI95upper, ymax=rg.CI95lower,width=.0), position=position_dodge(.7),stat="identity",size=1.3) +
  geom_point(position=position_dodge(.7), stat="identity",size=4,color="lightgrey") + #) +
  geom_text(aes(x=region,y=sig_val,label=sig),
            position=position_dodge(width=0.9),size=8, stat="identity") +
  
  # legend and axes
  labs(y=bquote(r[g]),shape=bquote(r[g] ~' Estimate type')) +
  scale_color_manual(values =colors,guide="none") +
  scale_fill_manual(values =colors,guide="none") +
  scale_shape_manual(values=c(24,21)) +
  scale_alpha_manual(values=c(0.5,1),guide="none") +
  scale_x_discrete(limits=rev(levels(rg_cogn$region))) +
  scale_y_continuous(breaks=seq(-0.6,0.6,by=0.30)) +
  # 
  coord_flip(ylim=c(-0.85,0.85)) +
  # coord_flip() +
  facet_grid( .~cognitive, switch="y",scales="free") +
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


rg_sibs_plot
leg_sibs<-get_legend(rg_sibs_plot + theme(legend.position = "right") +
                       guides(shape=guide_legend(override.aes = list(size=3,color="black"))) +
                       theme(legend.justification="center",legend.text=element_text(size=14),
                             legend.title=element_text(size=16))
                     )

rg_cogn_sibs_plot<-plot_grid(rg_cogn_plot2,
          plot_grid(rg_sibs_plot,leg_sibs,rel_widths = c(1,1), nrow=1),
          rel_heights = c(1,1),
          labels=c("a","b"),
          ncol=1)


brain_rgL_cogn_sibs_plot<-plot_grid(brains,rg_cogn_sibs_plot,nrow=2,rel_heights = c(1,2))
ggsave(brain_rgL_cogn_sibs_plot,file=paste0(out_dir,"brainL_rg_leftROIs_cognANDreadANDsibs_ukb_big40.png"),
       width=12,height=12)
