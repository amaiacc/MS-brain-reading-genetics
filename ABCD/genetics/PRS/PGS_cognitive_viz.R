# ** to add data etc, continuing after PGS_lmm.R **
rm(list=ls())
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects/resources/datasets/ABCD/scripts/genetics/PRS/base_cognitive.config"


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
PGS_dir=paste(dir,project,"/data/working_data/genotyping/prsice/","cognitive/",sep="")
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

#------------------------------------
# visualize
#------------------------------------
## COGN on BEHAVIOUR
max_p<- (-log10(subset(pgs_summary,subset==s)$P)) %>% max() %>% round(digits=0) +5
ea_PGS_plot_c<-ggplot(data=subset(pgs_summary,base_pheno=="EA3"&subset==s) %>%  filter(phenotype %in% vars2test_c),
                      aes(x=threshold,y=R2*100,fill=-log10(P))) +
  geom_col(width=1,position="dodge") + 
  # facet_grid(.~pheno,scale="free") +
  scale_fill_distiller(palette = 3, direction=1, name=expression(-log[10]*"(P-value)"),
                       limit=c(0,max_p),n.breaks=10,
                       guide = guide_legend(direction = "vertical", title.position = "right", 
                                            title.theme = element_text(angle = 90),title.vjust=0.5) ) +
  labs(y=expression("% variance " (Delta~R^2)),x="GWAS P-value threshold",
       title="PGS Educational Attainment") +
  # mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim=c(0,4)) +
  NULL


cp_PGS_plot_c<- ggplot(data=subset(pgs_summary,base_pheno=="CP"&subset==s) %>%  filter(phenotype %in% vars2test_c),
                      aes(x=threshold,y=R2*100,fill=-log10(P))) +
  geom_col(width=1,position="dodge") + 
  # facet_grid(.~target_pheno,scale="free") +
  scale_fill_distiller(palette = 3, direction=1, name=expression(-log[10]*"(P-value)"),  
                       limit=c(0,max_p),n.breaks=10,
                       guide = guide_legend(direction = "vertical", title.position = "right", 
                                            title.theme = element_text(angle = 90),title.vjust=0.5) ) +
  labs(y=expression("% variance " (Delta~R^2)),x="GWAS P-value threshold",
       title="PGS Cognitive Performance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim=c(0,4)) +
  NULL


thrlist=levels(pgs_summary$threshold)


## define best base dataset and p-value threshold for predicting reading
reading_best<-subset(pgs_summary,subset==s) %>%  filter(phenotype %in% vars2test_c)
reading_best<-reading_best %>% filter(R2==max(reading_best$R2))

best_base<-reading_best$base_pheno
best_thr<-reading_best$threshold

cp_PGS_plot_c2 <- cp_PGS_plot_c + geom_rect(aes(
  xmin=which(thrlist==best_thr)-0.5,
  xmax=which(thrlist==best_thr)+0.5,
  ymin =0, ymax = max(reading_best$R2)*100),
  color="black",size=1.5) +
  guides(fill=guide_legend(override.aes = list(color = "white"),
                           direction = "vertical", 
                           title.position = "right", 
                           title.theme = element_text(angle = 90),title.vjust=0.5)) +
  theme(axis.text.x=
          element_text(face=if_else(levels(pgs_summary$threshold) %in% as.character(best_thr),"bold","plain"),
                       size=if_else(levels(pgs_summary$threshold) %in% as.character(best_thr),11,10))
  ) +
  NULL

  
PGS_legend<-get_legend(ea_PGS_plot_c)

PGS_all_plot_c<-plot_grid(ea_PGS_plot_c + theme(legend.position = "none"),
                          cp_PGS_plot_c2 + theme(legend.position = "none"),
                          NULL,
                          PGS_legend,
                          nrow=1,
                          labels=c("A","B",""),rel_widths = c(1,1,0.05,0.2))


#---------------------------------------------------
# save summary plot and table
#---------------------------------------------------
ggsave(PGS_all_plot_c,file=paste0(out_dir,"PGS_effects_cognitive2ABCD_",s,"_behav.pdf"),width=11,height=3)

# clean intermediate
rm(list=ls(pattern="plot"))

#---------------------------------------------------

## COGN on BRAIN

pgs_brain<-subset(pgs_summary,subset==s) %>% filter(!phenotype %in% vars2test_c) %>%
  filter(base_pheno==best_base & threshold==best_thr) %>% distinct()

# remove adjusted if measure is total
w<-intersect(which(pgs_brain$adj!=""),grep("total",pgs_brain$target_pheno))
if(length(w)>0){
  pgs_brain<-pgs_brain[-w,]
}
rm(w)

#---------------------------------------------------
## individual level data
pgs_d<-readRDS(paste0(out_dir,"PGS_",base_project,"_",s,"_data.rds")) %>% distinct()

# plots of interest
# inspect specific relationships
plot_base_thr_pheno<-function(base_name,base_thr,p,cov="",fillcol="black"){
  PGS=paste(base_name,base_thr,sep="_")
  if(cov==""){
    tmp<-pgs_d[,c(p,PGS)] %>% distinct()
    colnames(tmp)<-c("pheno","PGS")
    # define shape for plot
    s=21
  } else {
    if(cov=="tCSA"){cov2="smri_area_cort.destrieux_total.lh"}
    tmp<-pgs_d[,c(p,PGS,cov2)] %>% distinct()
    tmp$res<-lm(paste(p,"~",cov2),data=tmp)$residuals
    tmp<- tmp[,c("res",PGS)]
    colnames(tmp)<-c("pheno","PGS")
    #define shape for plot
    s=24
    
  }
  
  tmp$pheno<-tmp$pheno %>% scale()
  tmp$PGS<-tmp$PGS %>% scale()

  p_2<-gsub("nihtbx_|_uncorrected|smri_area_cort.destrieux_","",p) %>% simpleCap()
  p_2<-if(cov!=""){paste(p_2,"\nadjusted for",cov)}

  # define  deciles for PGS, and plot
  tmp <- tmp %>% filter(!is.na(pheno)) %>% mutate(PGS_decile = ntile(PGS, 10)) %>%
    group_by(PGS_decile)
  r2<-subset(pgs_summary,phenotype==p&base_pheno==base_name&threshold==base_thr&adj==cov) %>% distinct() %>% pull(R2) %>% round(.,digits=4)
  z <- tmp %>%  summarise(
    mean = mean(pheno,na.rm=TRUE),
    sd = sd(pheno,na.rm=TRUE),
    se= (sd/sqrt(length(pheno)))
  )
  # define plot limits
  
  ymin=(min(z$mean) - 1.96*min(z$se)) %>% round(digits=1) - 0.2
  ymax=(max(z$mean) - 1.96*max(z$se)) %>% round(digits=1) + 0.2
  
  # plots
  p1 <- ggplot(data=z,aes(x=as.factor(PGS_decile),y=mean)) +
    # geom_dotplot(data=tmp,aes(x=as.factor(PGS_decile),y=pheno),binaxis='y',stackdir='center', 
    #              dotsize=0.3,alpha=0.2,binwidth=0.1,fill=fillcol) +
    geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), colour="black", width=.1) +
    geom_point(size=5, shape=s, fill=fillcol,alpha=alpha) +
    # geom_text(aes(x = 2, 
    #               y = ymax-0.05, #max(tmp$pheno,na.rm=TRUE), # max(mean+sd,na.rm=TRUE), 
    #               label = paste("Delta~R^2: ",r2 )) ,
    #           color = "black", size=3,
    #           parse = TRUE) +
    labs(y=p_2,x=paste("PGS-",base_name," ",base_thr," decile",sep="")) +
    ylim(c(ymin,ymax)) +
    NULL
  p2 <- ggplot(data=tmp,aes(x=as.factor(PGS_decile),y=pheno)) +
    geom_dotplot(binaxis='y',stackdir='center', dotsize=0.3,alpha=0.2,binwidth=0.1) +
    geom_text(aes(x = 1,
                  y = max(tmp$pheno,na.rm=TRUE),
                  label = paste("Delta~R^2: ",r2 )) ,
              color = "black",
              parse = TRUE) +
    labs(y=p_2,x=paste("PGS-",base_name," ",base_thr," decile",sep="")) +
    NULL
  p3 <- ggplot(data=tmp,aes(x=PGS,y=pheno)) + geom_point(alpha=0.2) +
    geom_smooth( method="lm",  se=TRUE, size=2) +
    geom_text(aes(x = min(tmp$PGS,na.rm=TRUE), 
                  y = max(tmp$pheno,na.rm=TRUE), 
                  label = paste("Delta~R^2: ",r2 )) ,
              color = "black",
              parse = TRUE) +
    labs(y=p_2,x=paste("PGS-",base_name," ",base_thr,sep="")) +
    NULL
  
  return(plot_grid(p1,p3,ncol=1))
}


plotdecile_base_thr_pheno<-function(base_name,base_thr,p,cov="",fillcol="black",alpha=1){
  PGS=paste(base_name,base_thr,sep="_")
  if(cov==""){
    tmp<-pgs_d[,c(p,PGS)] %>% distinct()
    colnames(tmp)<-c("pheno","PGS")
    # define shape for plot
    s=21
  } else {
    if(cov=="tCSA"){cov2="smri_area_cort.destrieux_total.lh"}
    tmp<-pgs_d[,c(p,PGS,cov2)] %>% distinct()
    tmp$res<-lm(paste(p,"~",cov2),data=tmp)$residuals
    tmp<- tmp[,c("res",PGS)]
    colnames(tmp)<-c("pheno","PGS")
    #define shape for plot
    s=24
    
  }
  
  # tmp$pheno<-tmp$pheno %>% scale()
  tmp$PGS<-tmp$PGS %>% scale()
  
  p_2<-gsub("nihtbx_|_uncorrected|smri_area_cort.destrieux_","",p) %>% gsub("\\.lh","",.) %>% gsub("\\."," ",.) %>% simpleCap()
  p_2<-if(cov!=""){paste(p_2,"\nadjusted for",cov)}else{p_2<-p_2}
  # define  deciles for PGS, and plot
  tmp <- tmp %>% filter(!is.na(pheno)) %>% mutate(PGS_decile = ntile(PGS, 10)) %>%
    group_by(PGS_decile)
  r2<-subset(pgs_summary,phenotype==p&
               base_pheno==base_name&
               threshold==base_thr&
               adj==cov) %>% distinct() %>% pull(R2) %>% round(.,digits=4)
  z <- tmp %>%  summarise(
    mean = mean(pheno,na.rm=TRUE),
    sd = sd(pheno,na.rm=TRUE),
    se= (sd/sqrt(length(pheno)))
  )
  
  # define plot limits
  
  ymin=(min(z$mean) - 1.96*min(z$se)) %>% round(digits=1) - 0.2
  ymax=(max(z$mean) - 1.96*max(z$se)) %>% round(digits=1) + 0.2
  
  
  p1 <- ggplot(data=z,aes(x=as.factor(PGS_decile),y=mean)) +
    # geom_dotplot(data=tmp,aes(x=as.factor(PGS_decile),y=pheno),binaxis='y',stackdir='center',
    #              dotsize=0.3,alpha=0.2,binwidth=0.1,fill=fillcol) +
    geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), colour="black", width=.1) +
    geom_point(size=5, shape=s, fill=fillcol,alpha=alpha) +
    geom_text(aes(x = 2, 
                  y = ymax-0.05, #max(tmp$pheno,na.rm=TRUE), # max(mean+sd,na.rm=TRUE), 
                  label = paste("Delta~R^2: ",r2 )) ,
              color = "black", size=3,
              parse = TRUE) +
    labs(y=p_2,x=paste("PGS-",base_name," ",base_thr," decile",sep="")) +
    ylim(c(ymin,ymax)) +
    NULL

  
  return(p1)
}

# for reading
# reading_best_plot<-plot_base_thr_pheno(base_name=best_base,base_thr=best_thr,p="nihtbx_reading_uncorrected")

reading_best_decileplot<-plotdecile_base_thr_pheno(base_name=best_base,base_thr=best_thr,p="nihtbx_reading_uncorrected")


# cross-trait
csa_best_decileplot<-plotdecile_base_thr_pheno(base_name=best_base,base_thr=best_thr,
                                               p="smri_area_cort.destrieux_total.lh",
                                               fillcol=braincolor_codes %>% filter(measure=="AREA"&Region=="total") %>% pull(color)) +
  ylab("Total CSA") +
  xlab(expression('PGS'[CP]))
gtempsup_best_decileplot<-plotdecile_base_thr_pheno(base_name=best_base,base_thr=best_thr,
                                                    p="smri_area_cort.destrieux_g.temp.sup.lateral.lh",
                                                    fillcol=braincolor_codes %>% filter(measure=="AREA"&Region=="g.temp.sup.lateral") %>% pull(color)) +
  ylab("STG Lateral CSA") +
  xlab(expression('PGS'[CP]))

gtempsupAdj_best_decileplot<-plotdecile_base_thr_pheno(base_name=best_base,base_thr=best_thr,
                                                    cov="tCSA",
                                                    p="smri_area_cort.destrieux_g.temp.sup.lateral.lh",
                                                    fillcol=braincolor_codes %>% filter(measure=="AREA"&Region=="g.temp.sup.lateral") %>% pull(color),
                                                    alpha=0.5) +
  ylab("STG Lateral CSA\nadjusted for Total CSA") +
  xlab(expression('PGS'[CP]))


CSA_all_decileplot<-plot_grid(gg_area_all2,csa_best_decileplot,
          gtempsup_best_decileplot,
          gtempsupAdj_best_decileplot)


CSA_all_decileplot2<-plot_grid(csa_best_decileplot,
                               gtempsup_best_decileplot,
                               gtempsupAdj_best_decileplot, nrow=1)


ggsave(reading_best_decileplot,file=paste0(out_dir,"PGS_CP0.05_decileplot_reading.pdf"),height=4,width=4)


ggsave(CSA_all_decileplot,file=paste0(out_dir,"PGS_CP0.05_decileplot_CSA_brain.pdf"),height=5,width=7)
ggsave(CSA_all_decileplot2,file=paste0(out_dir,"PGS_CP0.05_decileplot_CSA.pdf"),height=3,width=12)

# out_dir2 to combine figures for manuscript
out_dir2=paste0(dir,"ABCD/data/working_data/output/")
saveRDS(CSA_all_decileplot2,file=paste0(out_dir2,"CSA_all_decileplot2.rds"))
