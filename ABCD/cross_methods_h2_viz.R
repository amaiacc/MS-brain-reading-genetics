# clean workspace
rm(list=ls())
# libraries and custom functions
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(RColorBrewer)
library(tidyverse)
#---------------------------------------------------
options(stringsAsFactors = FALSE)
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
# read ggseg images and color table, for consistency
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/",sep="")
load(paste(pheno_dir,"ggseg_brain_ROIs.Rdata",sep=""))
braincolor_codes<- braincolor_codes %>% 
  mutate(region_name=Region  %>% gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% toupper() )
# rois colors
colors_a <- braincolor_codes  %>% filter(hemisphere=="L" & measure=="AREA") %>% pull(color)
colors_t <- braincolor_codes  %>% filter(hemisphere=="L" & measure=="THICKNESS") %>% pull(color)
colors<- braincolor_codes %>% pull(color) %>% unique()  
#----------------------------------------------------------------------
# set working dir
working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/",sep="")
setwd(working_dir)
#
out_dir=paste0(working_dir,"/output/")
if(!dir.exists(out_dir)){dir.create(out_dir)}


# define data dirs
gcta_dir=paste0(working_dir,"/genotyping/GCTA/output/")
ldsc_dir=paste0(working_dir,"/genotyping/ldsc/")
ldsc_dir2=paste0(dir,"/resources/datasets/GWAS_sumstats/ldsc/")
ldsc_dir3=paste0(dir,"coeduca/data/working_data/GenLang/projects/GWASMA_freeze_v1/data/ldsc/",sep="")
#----------------------------------------------------------------------
# read data
## ldsc, from GWAS sumstats
# brain ldsc
ldsc1_h2<-read.csv(paste0(ldsc_dir,"ldsc_h2_UKB_BIG40.csv")) %>% 
  mutate(region_name=factor(region,levels=unique(braincolor_codes$region_name))) %>%
  mutate(p1=IDP_short_name,
         pval=p,
         dataset="GWAS sumstats",
         method="LDSC",
         category="Brain") %>%
  mutate(adj=if_else(region_name=="TOTAL","-","adjGlobal")) %>%
  dplyr::select(p1,h2,h2.CI95upper,h2.CI95lower,pval,category,hemisphere,measure,region_name,dataset,method,adj)
# cognitive ldsc
ldsc2_h2<-read.csv(paste0(ldsc_dir2,"ldsc_h2_publicGWASes.csv")) %>%
  filter(Trait_name=="EA"|Trait_name=="CP"|Trait_name=="IQ"|Trait_name=="DDyslexia") %>%
  mutate(p1=Trait_name,
         pval=p,
         dataset="GWAS sumstats",
         method="LDSC",
         category="Cognitive") %>%
  dplyr::select(p1,h2,h2.CI95upper,h2.CI95lower,pval,category,dataset,method)
# combine all ldsc
ldsc_h2<-merge(ldsc1_h2,ldsc2_h2,all=TRUE)
# ldsc_h2<-merge(ldsc_h2,ldsc3_h2,all=TRUE)
rm(ldsc1_h2,ldsc2_h2)
# rm(ldsc3_h2)
## gcta
gcta_h2<-read.csv(paste0(gcta_dir,"gcta_estimates_h2_ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned_rm05_adj.csv")) %>% 
  mutate(region_name=region %>% gsub("smri_thick_cort.|smri_area_cort.|destrieux_","",.) %>% gsub("_adjGlobal","",.) %>%
           gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% gsub("AREA_|area_|THICK_|thick_","",.) %>% toupper()) %>%
  mutate(h2.CI95upper=h2+1.96*se,
         h2.CI95lower=h2-1.96*se,
         dataset="ABCD",
         method="GREML") %>%
  dplyr::select(p1,h2,h2.CI95upper,h2.CI95lower,pval,category,hemisphere,measure,region_name,dataset,method,adj)

## combine all
h2<-merge(gcta_h2,ldsc_h2,all=TRUE)
rm(gcta_h2,ldsc_h2)

h2$method<-factor(h2$method,levels=c("GREML","LDSC"))
h2$p<-as.numeric(h2$pval)
h2$p[h2$pval=="< 0.001"]<-0.001

#
h2 <- h2 %>% mutate(name=if_else(str_detect(p1,"WR|reading"),"Word Reading",
                    if_else(str_detect(p1,"IQ|fluid"),"Fluid",
                        if_else(str_detect(p1,"wiscv"),"Matrix reasoning",
                            if_else(str_detect(p1,"CP"),"Cognitive Performance",
                                    if_else(str_detect(p1,"cryst"),"Crystalized",
                                            if_else(str_detect(p1,"picvocab"),"Vocabulary",
                                                    if_else(str_detect(p1,"Dyslexia"),"Dyslexia",
                                                            if_else(str_detect(p1,"EA"),"Educational Attainement","-"
                                                            ))))))))) %>%
  mutate(name2=if_else(str_detect(p1,"WR|reading|Dyslexia"),"Reading",
                       if_else(str_detect(p1,"IQ|CP|cryst|fluid|wisc"),"Intelligence",
                               if_else(str_detect(p1,"picvocab"),"Vocabulary",
                                       if_else(str_detect(p1,"EA"),"Education","-"
                                       ))))) 
w<-which(is.na(h2$measure))
h2$measure[w]<-h2$name2[w]
h2$region_name[w]<-h2$name[w]
rm(w)
# subset and clean columns
h2<-subset(h2,measure!="-"&(hemisphere=="L"|is.na(hemisphere)))
h2<-h2 %>% mutate(
  hemisphere=if_else(is.na(hemisphere)|hemisphere=="-","",hemisphere),
  adj=if_else(is.na(adj)|adj=="-","",adj),
  measure=if_else(measure=="AREA","CSA",
                  if_else(measure=="THICKNESS","CT",measure)))

# save table with all estimates
h2 %>% 
  mutate(trait=
           paste(region_name," ",hemisphere," (",paste(measure,adj,sep=" "),")",sep="") %>% 
           gsub("  ", " ",. ) %>% gsub(" $", "",. ) %>% gsub(" )", ")",. ) 
         ) %>%
  dplyr::select(trait,dataset,method,category,measure,region_name,hemisphere,adj,
         h2,h2.CI95upper,h2.CI95lower,pval) %>%
  arrange(-h2) %>%
  rename(adjustement=adj, name=region_name) %>% 
  write.csv(.,file=paste0(out_dir,"h2_cross_methods.csv"),row.names = FALSE)

#----------------------------------------------------------------------
# plots
#----------------------------------------------------------------------
## def colors
mytheme<-  theme(axis.title.x=element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1,face="bold"),
                 axis.text=element_text(size=12),
                 legend.position="bottom")
## BRAIN ##
h2_brain<-subset(h2,category!="Cognitive") %>% filter(hemisphere=="L") %>%
  mutate(region_name=factor(region_name,levels=unique(braincolor_codes$region_name))) %>% arrange(region_name) %>%
  mutate(region_name2=factor(region_name %>% as.character() %>% tolower() %>% sapply(.,simpleCap),
                             levels=unique(braincolor_codes$region_name)%>% tolower() %>% sapply(.,simpleCap) )
         )
h2_brain$measure<-h2_brain$measure %>% tolower() %>% sapply(.,simpleCap)
# add sig flag for viz
h2_brain$sig<-""
h2_brain$sig[h2_brain$p<0.05]<-"*"
h2_brain$sig[h2_brain$p<(0.05/(length(unique(h2_brain$region_name))-1))]<-"***"

## dummy var to code shapes similar to other plots (were method is not a variable)
h2_brain$method<-factor(h2_brain$method,levels=c("ACE","GREML","LDSC"))
h2_brain$shape<-paste(h2_brain$method,h2_brain$adj,sep="_") %>% 
  factor(., levels=paste(rep(c("ACE","GREML","LDSC"),each=2),c("-","adjGlobal"),sep="_"))
##
h2_brain_plot <- ggplot(data=h2_brain,
                        aes(x=region_name2,y=h2,width=0.7,
                            color=region_name,#fill=region_name,
                            shape=method
                            # shape=shape
                            )) +
  # annotate area and thickness, and add line in between
  annotate("text", x=7-0.25,y=1,label="CSA",size=3) +
  annotate("text", x=6+0.25,y=1,label="CT",size=3) +
  geom_vline(xintercept = seq(1.5,9.5,by=1), color="grey",linetype="solid") +
  geom_vline(xintercept = 7-0.5,color="black",alpha=0.8,linetype="solid",size=1.2) +
  #
  geom_hline(yintercept = 0,color="black",alpha=0.3,linetype="dashed") +
  geom_hline(yintercept = 1,color="black",alpha=0.3,linetype="dashed") +
  # geom_bar(aes(y=h2,alpha=hemisphere),position=position_dodge(), stat="identity",size=1) +
  geom_errorbar(aes(ymin=h2.CI95upper , ymax=h2.CI95lower,width=0), position=position_dodge(.7),size=1,alpha=0.7) +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_text(aes(x=region_name2,y=h2.CI95upper+0.2,label=sig),
            position=position_dodge(width=0.7),size=6, stat="identity",color="black",alpha=0.7) +
  ylab(bquote(''*h^2*'')) +
  scale_color_manual(values =colors,guide="none") +
  scale_fill_manual(values =colors,guide="none") +
  scale_shape_manual(name="Method",values=c(21,8,18)) +
  # scale_shape_manual(name="Method",values=c(19,17,1,2,13,4)) +
  scale_y_continuous(limits=c(-0.10,1.2), breaks = seq(0, 1, by = 0.20)) +
  # # scale_x_discrete(labels=h2_brain$region_name2 %>% unique()) +
  scale_x_discrete(limits=rev(levels(h2_brain$region_name2))) +
  # theme(axis.title.x=element_blank(),legend.position="bottom") +
  # mytheme +
  # NULL
  facet_grid(. ~adj, scales="free",space="free") +
  coord_flip() +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  theme(axis.title.y=element_blank(),legend.position="bottom") +
  NULL

## COGNITIVE ##
h2_cogn <- subset(h2,category=="Cognitive") %>% distinct() %>%
  filter(name!="-") %>% filter(p1!="nihtbx_cryst_uncorrected") %>%
  mutate(name2=factor(name2,levels=c("Reading","Vocabulary","Intelligence","Education"))) %>% 
  mutate(name=factor(name,levels=c("Word Reading","Dyslexia","Vocabulary",
                                   "Crystalized",
                                   "Fluid","Matrix reasoning",
                                   "Cognitive Performance",
                                   "Educational Attainement"))) %>% 
  arrange(name2,name)

# add sig flag for viz
h2_cogn$sig<-""
h2_cogn$sig[h2_cogn$p<0.05]<-"*"
h2_cogn$sig[h2_cogn$p<0.05/length(unique(h2_cogn$name))]<-"***"

# 
h2_cogn_plot <- ggplot(data=h2_cogn,aes(x=name,y=h2,width=0.7,color=name2,shape=method)) +
  # geom_vline(xintercept = c(2.5), color="grey",linetype="solid") +
  geom_hline(yintercept = c(0,1),color="black",alpha=0.3,linetype="dashed") +
  geom_errorbar(aes(ymin=h2.CI95upper , ymax=h2.CI95lower,width=0), position=position_dodge(.7),size=1,alpha=0.7) +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_text(aes(x=name,y=h2.CI95upper+0.2,label=sig),position=position_dodge(width=0.7),size=6,color="black", stat="identity",alpha=0.7) +
  ylab(bquote(''*h^2*'')) +
  scale_color_brewer(palette="Set2",guide="none") +
  scale_fill_brewer(palette="Set2",guide="none") +
  scale_shape_manual(name="Method",values=c(21,8,18)) +
  # facet_grid(.~ name2,scales = "free_x",space="free_x") +
  facet_grid( name2  ~., switch="y",scales="free",space="free") +
  scale_y_continuous(limits=c(-0.10,1.2), breaks = seq(0, 1, by = 0.20)) +
  scale_x_discrete(drop = TRUE) +
  coord_flip() +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  theme(axis.title.y=element_blank(),legend.position="bottom") +
  NULL



legend_h2 <- get_legend(  h2_cogn_plot +
                           # guides(shape = guide_legend(nrow = 1)) +
                            theme(legend.position = "right",legend.direction="vertical") )

# now add the title
title <- ggdraw() + 
  draw_label(
    "Heritability",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

h2_plot<-plot_grid(
  plot_grid(h2_cogn_plot + theme(legend.position="none") + labs(title="Cognitive measures"),
            NULL,
            legend_h2,
            NULL,
            nrow=1,
            rel_widths=c(2,0.2,0.3,0.5)),
  # h2_cogn_plot + theme(legend.position="right") + labs(title="Cognitive measures"),
  h2_brain_plot + theme(legend.position="none") + labs(title="Reading associated brain measures"),
  rel_heights = c(1,1.2),
  ncol=1,  #align="v",
  labels=c("A","B"),
  axis ="l"
)

h2_plot2<-plot_grid(
  plot_grid(h2_cogn_plot + theme(legend.position="none") + labs(title="Cognitive measures"),
            NULL,
            legend_h2,
            NULL,
            nrow=1,
            rel_widths=c(2,0.2,0.3,0.5)),
  # h2_cogn_plot + theme(legend.position="right") + labs(title="Cognitive measures"),
  h2_brain_plot + theme(legend.position="none") + labs(title="Reading associated brain measures"),
  rel_heights = c(1,1.2),
  ncol=1,  #align="v",
  labels=c("",""),
  axis ="l"
)
# save
ggsave(h2_brain_plot,file=paste0(out_dir,"h2_brain_cross_methods.png"))
ggsave(h2_cogn_plot,file=paste0(out_dir,"h2_cognitive_cross_methods.png"))
ggsave(h2_plot,file=paste0(out_dir,"h2_all_cross_methods.pdf"),width=8,height=11)
ggsave(h2_plot2,file=paste0(out_dir,"h2_all_cross_methods_nolabels.pdf"),width=8,height=11)

ggsave(h2_plot,file=paste0(out_dir,"h2_all_cross_methods.png"),width=8,height=11)
ggsave(h2_plot2,file=paste0(out_dir,"h2_all_cross_methods_nolabels.png"),width=8,height=11)



# out_dir2 to combine figures for manuscript
out_dir2=paste0(dir,"resources/datasets/ABCD/data/working_data/output/")
saveRDS(h2_plot2,file=paste0(out_dir2,"h2_all_cross_methods_nrow1.rds"))
saveRDS(h2_plot,file=paste0(out_dir2,"h2_all_cross_methods_ncol1.rds"))

