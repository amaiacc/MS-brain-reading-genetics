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
working_dir=paste(base_dir,"rois2followup_RH_","baseline","/",sep="")
out_dir=paste(base_dir,"/summary/figures/",sep="")


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
rois2follow <-read.table("summary/tables/ROImeasures2follow_RH.txt") %>% 
   mutate(model=gsub("lh","rh",V1)) %>% select(model)

# baseline models
ttables<-c( paste(paste(working_dir,"/tables/",sep=""),
                  list.files(paste(working_dir,"/tables/",sep=""),pattern="t_table"),sep="/"))

ttables<-ttables[grep("m0b",ttables,invert=TRUE)] # not TBV corrected

for (f in ttables){
  
  out_dir=strsplit(f,"tables") %>% sapply("[[",1) %>% paste(.,"/figures/",sep="")
  if(!dir.exists(out_dir)){dir.create(out_dir)}
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
  
  # make figures
  for (m in c(unique(t$measure)) ){
    
    titlep <- paste(dv, " ~ covariates ", c, " + ", m, " ROI",sep="") %>% gsub("  "," ",.) # "(",hemi,")"
    t2 <- t %>% filter(measure==m &hemi==hemi&region!="total"&region!="mean") %>% filter(P<0.05)
    
    p <- t2 %>% select(label,t.value,P,model) %>% brain_join(desterieux) %>% filter(hemi=="right") %>%
      ggplot() + 
      geom_sf(aes(fill=t.value),colour="black",size=0.5) +
      scale_fill_gradient2(name="T-value",limit=c(-10,10)) +
      theme_void() +
      labs(caption=titlep) +
      NULL

    assign(paste("plot_",m,"_m",i,sep=""),p)
    rm(p)
    rm(t2)
  }
  rm(m,titlep)
}
rm(t,i,c,f)

#
d$region[is.na(d$region)]<-d$region.name[is.na(d$region)]
#---------------------
# plots
# make figures combining different models
thick_list<-ls(pattern="plot_thick")
area_list<-ls(pattern="plot_area")
# extract a legend that is laid out horizontally
legend_t <- get_legend(  get(thick_list[1]) +
                           guides(color = guide_legend(nrow = 1)) +  theme(legend.position = "bottom") )
legend_a <- get_legend(  get(area_list[1]) +
                           guides(color = guide_legend(nrow = 1)) +  theme(legend.position = "bottom") )
#
n_t<-ceiling(length(thick_list[grep("sig",thick_list,invert=TRUE)])/4)
n_a<-ceiling(length(area_list[grep("sig",thick_list,invert=TRUE)])/4)


#---------------------
# visualize estimates, across runs
library(RColorBrewer)

# get order of estimates, from the non adjusted LM:
betas_order_csa <- d %>% filter(adj==""&measure=="area") %>% arrange(Estimate)
betas_order_ct <- d %>% filter(adj==""&measure=="thick") %>% arrange(Estimate)

d_csa <- d %>% filter(measure=="area") %>% filter(!(region=="total"&adj!="")) %>% 
  mutate(region=factor(region,levels=betas_order_csa %>% pull(region) %>% unique())) %>% 
  mutate(model=factor(model,levels=betas_order_csa %>% pull(model) %>% unique()))

d_ct<-d %>% filter(measure=="thick")%>% filter(!(region=="mean"&adj!="")) %>% 
  mutate(region=factor(region,levels=betas_order_ct %>% pull(region) %>% unique())) %>% 
  mutate(model=factor(model,levels=betas_order_ct %>% pull(model) %>% unique()))

csa_betas <- ggplot(data=d_csa,aes(x=Estimate,y=region,color=adj,shape=adj)) + 
  geom_vline(xintercept = 0,color="black",alpha=0.7,linetype="solid") +
  # highlight specific regions that will be followed up
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = "S temporal sup", ymax = Inf),
  #           color="white",fill = "lightgrey", alpha = 0.05) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = "S occipital ant", ymax = "S precentral-sup-part"),
  #           color="white",fill = "lightgrey", alpha = 0.05) +
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
  labs(title="",y="",color="Model",shape="Model") +
  # same x range for csa and ct
  coord_cartesian(xlim=c(min(d$Estimate)-0.02,max(d$Estimate)+0.02) %>% round(digits=2)) +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  theme(axis.title.y=element_blank(),legend.position="bottom",legend.direction = "vertical") +
  NULL

ct_betas <- ggplot(data=d_ct,aes(x=Estimate,y=region,color=adj,shape=adj)) + 
  geom_vline(xintercept = 0,color="black",alpha=0.7,linetype="solid") +
  # highlight specific regions that will be followed up
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = "S calcarine", ymax = Inf),
  #           color="white",fill = "grey", alpha = 0.05) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = "S subparietal"),
  #           color="white",fill = "grey", alpha = 0.05) +
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
  labs(title="",y="",color="Model",shape="Model") +
  # same x range for csa and ct
  coord_cartesian(xlim=c(min(d$Estimate)-0.02,max(d$Estimate)+0.02) %>% round(digits=2)) +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  theme(axis.title.y=element_blank(),legend.position="bottom",legend.direction = "vertical") +
  NULL
  

all_betas<-plot_grid(
  plot_grid(NULL,csa_betas,ncol=1,rel_heights = c(2,3.5)),
  ct_betas,
  nrow=1, labels=c("A","B"))

## betas
ggsave(csa_betas,file=paste(out_dir,pheno,"_",atlas,"_","area_betas",".pdf",sep=""),height=4,width=7)
ggsave(ct_betas,file=paste(out_dir,pheno,"_",atlas,"_","thick_betas",".pdf",sep=""),height=4,width=7)
ggsave(all_betas,file=paste(out_dir,pheno,"_",atlas,"_","all_betas",".pdf",sep=""),height=6,width=15)

#---------------------
Ttitle <- ggdraw() + draw_label("CT (RH)", fontface='bold')
Atitle <- ggdraw() + draw_label("CSA (RH)", fontface='bold')

plot_thick_m0d_all <- plot_grid( Ttitle, plot_grid(plotlist =lapply(thick_list[grep("m0d$|m0d_ROIs",thick_list,invert=FALSE)], 
                                                          function(x) {get(x) + theme(legend.position="none")+labs(caption="")}),
                                         ncol=1),  legend_t , ncol = 1, rel_heights = c(.1, 1, .1)) 
plot_area_m0c_all <- plot_grid( Atitle, plot_grid(plotlist =lapply(area_list[grep("m0c$|m0c_ROIs",area_list,invert=FALSE)], 
                                                          function(x) {get(x) + theme(legend.position="none")+labs(caption="")}),
                                         ncol=1),  legend_t , ncol = 1, rel_heights = c(.1, 1, .1)) 



# combined plots
out_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/","/output/")

if(!dir.exists(out_dir)){dir.create(out_dir)}

figure <- plot_grid(
          plot_grid(NULL,csa_betas,ncol=1,rel_heights = c(2,3.5),align="hv"),
          plot_grid(
            plot_area_m0c_all +  theme(plot.background = element_rect(color = brewer.pal(12, "Paired")[3],size=2)),
            plot_thick_m0d_all +  theme(plot.background = element_rect(color = brewer.pal(12, "Paired")[7],size=2)),
            ncol=1,rel_heights = c(1,1),
            labels=c("B","C")),
          ct_betas,
          nrow=1,
          # rel_widths = c(1,1,1),
          labels=c("A","","D")
          )


ggsave(figure,file=paste(out_dir,pheno,"_",atlas,"_rois2follow_RH",".pdf",sep=""),width=20,height=6)

#---------------------
# clean intermediate
rm(list=thick_list)
rm(list=area_list)
rm(thick_list,area_list,ttables)
rm(n_t,legend_t,n_a,legend_a)

