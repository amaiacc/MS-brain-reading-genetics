# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/reading.config" # parameters for this run
# config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/picvocab.config" # parameters for this run
#---------------------------------------------------
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library("ggseg")
library("ggsegExtra")
# ggseg::install_atlases()
library(ggsegDesterieux)
library(ggsegDefaultExtra)
#------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
# source config file to define parameters for this run
source(config_file)
source(paste0(scripts_dir,"color_functions.R"))

#---------------------
if(atlas=="destrieux"){
  target_labels<-data.frame(target=ggseg(atlas=desterieux)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
} else if (atlas=="desikan"){
  target_labels<-data.frame(target=ggseg(atlas=dkt)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
}

base_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",atlas,"/",sep="")
setwd(base_dir)
# read ROIs to follow from brain-behaviour associations
summary_codes<-read.table("summary/tables/ROImeasures2follow.txt") %>% 
  rename(model=V1) %>%
  mutate(region_h=strsplit(model,atlas) %>% sapply("[[",2) %>% gsub("^_","",.),
         region= region_h %>% gsub(".lh|.rh","",.),
         hemisphere=if_else(str_detect(model,"lh"),"L",if_else(str_detect(model,"rh"),"R","NA")),
         measure=if_else(str_detect(model,"area"),"AREA",if_else(str_detect(model,"thick"),"THICKNESS","NA"))
         )


# remove unnecessary variables from config file, to keep the working space clean
rm(args,cov11,cov12,cov2,covs_int,thr, trim_val,gen_v,dv,pc_file)
rm(s,v)

#------------------------------------------
# plot
braincolor_codes<-colors_brain(data=summary_codes,atlas=atlas,p_col="region_h") 
## make ggseg plots for color coding regions
titleA <- ggdraw() + 
  draw_label("CSA",fontface = 'bold',x = 0,hjust = 0 ) +   theme(    plot.margin = margin(0, 0, 0, 7) )
titleT <- ggdraw() + 
  draw_label("CT",fontface = 'bold',x = 0,hjust = 0 ) +   theme(    plot.margin = margin(0, 0, 0, 7) )

gg_ROIarea <- braincolor_codes %>% filter(measure=="AREA"&!is.na(label)) %>%  
  select(label,measure,color,Region) %>% brain_join(desterieux) %>%
  filter(hemi=="left") %>%
  ggplot() +
  geom_sf(aes(fill=Region),color="lightgrey",size=0.5) +
  theme_void() +
  scale_fill_manual(values=braincolor_codes %>% filter(measure=="AREA"&!is.na(label)&hemi=="left") %>% pull(color),
                    labels=c(braincolor_codes %>% filter(measure=="AREA"&!is.na(label)&hemi=="left") %>% pull(Region) %>% as.character(),"n.s.")
                    ) +
  theme(legend.position="none") +
  NULL
  
gg_ROIarea2 <- braincolor_codes %>% filter(measure=="AREA"&!is.na(label)) %>% filter(Region=="g.temp.sup.lateral") %>%  
  select(label,measure,color,Region) %>%  brain_join(desterieux) %>%
  filter(hemi=="left") %>%
  ggplot() +
  geom_sf(aes(fill=Region),color="lightgrey",size=0.5) +
  theme_void() +
  scale_fill_manual(values=braincolor_codes %>% filter(measure=="AREA"&!is.na(label)&hemi=="left") %>% filter(Region=="g.temp.sup.lateral") %>% pull(color),
                    labels=c(braincolor_codes %>% filter(measure=="AREA"&!is.na(label)&hemi=="left") %>% filter(Region=="g.temp.sup.lateral") %>% pull(Region) %>% as.character(),"n.s.")) +
  theme(legend.position="none") +
  NULL

gg_area <-  braincolor_codes %>% filter(measure=="AREA"&!is.na(label)) %>%  
  select(label,measure,color,Region) %>%  brain_join(desterieux) %>%
  filter(hemi=="left") %>%
  ggplot() +
  geom_sf(fill=braincolor_codes %>% filter(measure=="AREA"&is.na(label)&hemi=="left") %>% pull(color),color="lightgrey",size=0.5) +
  # scale_fill_manual(values=braincolor_codes %>% filter(measure=="AREA"&is.na(label)&hemi=="left") %>% pull(color), labels="Total") +
  theme_void() +
  theme(legend.position="none") +
  NULL

gg_area_all<-plot_grid(gg_area,gg_ROIarea,ncol=1)

gg_area_all2<-plot_grid(gg_area,gg_ROIarea2,ncol=1)


gg_ROIthick <- braincolor_codes %>% filter(measure=="THICKNESS"&!is.na(label)) %>%  
  select(label,measure,color,Region) %>% brain_join(desterieux) %>%
  filter(hemi=="left") %>%
  ggplot() +
  geom_sf(aes(fill=Region),color="lightgrey",size=0.5) +
  theme_void() +
  scale_fill_manual(values=braincolor_codes %>% filter(measure=="THICKNESS"&!is.na(label)&hemi=="left") %>% pull(color),
                    labels=c(braincolor_codes %>% filter(measure=="THICKNESS"&!is.na(label)&hemi=="left") %>% pull(Region) %>% as.character(),"n.s.")
  ) +
  theme(legend.position="none") +
  NULL


gg_thick_all<-plot_grid(NULL,gg_ROIthick,ncol=1)

# extract legends
# legend1<-get_legend(gg_area + guides(fill=guide_legend(title="Total Area")) + theme(legend.position="bottom",legend.direction="horizontal") )
legend2<-get_legend(gg_ROIarea + theme(legend.position="bottom") + 
                      guides(fill=guide_legend(title="CSA",ncol=3,direction="vertical"))+
                      theme(legend.justification="center")
)
legend3<-get_legend(gg_ROIthick + theme(legend.position="bottom") +  
                      guides(fill=guide_legend(title="CT",ncol=3,direction="vertical")) +
                      theme(legend.justification="center")
                    )
legend<-plot_grid(legend2,legend3,ncol=2,align="hv")
legend<-plot_grid(legend2,legend3,ncol=1,align="v")


gg_legend_thick_all<-plot_grid(legend,gg_ROIthick,ncol=1)


#----------------------------------------------------------------------
# remove intermediate objects
rm(target_labels)
rm(dir,base_dir,config_file)
#----------------------------------------------------------------------
# save plots, table and object
setwd("summary")
# 
write.csv(summary_codes,"ggseg_brain_colors4ROIs.csv",row.names = FALSE)

## brains
ggsave(gg_area,file="gg_areaTotal.png",height = 3,width = 7)
ggsave(gg_ROIarea,file="gg_ROIarea.png",height = 3,width = 7)
ggsave(gg_ROIthick,file="gg_ROIthick.png",height = 3,width = 7)
ggsave(gg_thick_all,file="gg_thick.png")
ggsave(gg_area_all,file="gg_area.png")
ggsave(gg_area_all2,file="gg_area_gen.png")

ggsave(plot_grid(
  plot_grid(titleA,gg_area_all,ncol=1,rel_heights = c(0.1,1),align="v"),
  plot_grid(titleT,gg_thick_all,ncol=1,rel_heights = c(0.1,1)),
  nrow=1,rel_widths = c(1,1)),
  file="gg_AreaThick.png",height=3, width=7)

# object
save.image(file="ggseg_brain_ROIs.Rdata")
