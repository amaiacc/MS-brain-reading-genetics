# clean workspace
rm(list=ls()); gc()
# define parameters from config file
config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/reading.config" # parameters for this run
hemi="lh"
#---------------------------------------------------
# define libraries
library(dplyr); library(tidyr)
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(ggrepel)
library(RColorBrewer)

# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
geno_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/")
scripts_dir=paste0(dir,"/resources/datasets/ABCD/scripts/phenotypes/")
# source config file to define parameters for this run
# source(paste(scripts_dir,config_file,sep="")) # if hard-coded
source(config_file) # as parameter

#
input_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/")
working_dir=paste0(input_dir,"/comparison","/",pheno,"_",atlas,"/")
if(!dir.exists(paste0(input_dir,"/comparison"))){dir.create(paste0(input_dir,"/comparison"))}
if(!dir.exists(working_dir)){dir.create(working_dir)}

# general functions to be used
source(paste(scripts_dir,"general_functions.R",sep=""))

# get to working dir
setwd(working_dir)
#---------------------------------------------------
# get subset definitions, as defined in define_subsets.R
## and subset to selected individuals only (i.e. s)
# id_subset<-read.table(paste0(input_dir,gen_v,"_subsets.table"),header=TRUE)
# # ids per subset
# ids_unrelated<-read.table(file=paste(input_dir,"/unrelated/","baseline_",pheno,"_",atlas,"/ABCDv2.0.1_update_DEAP_baseline_sMRI_baseline_subjectIDs_","unrelated","_trim","0.02",".txt",sep="")) 
# ids_eur<-read.table(file=paste(input_dir,"/eur/","baseline_",pheno,"_",atlas,"/ABCDv2.0.1_update_DEAP_baseline_sMRI_baseline_subjectIDs_","eur","_trim","0.02",".txt",sep="")) 
#---------------------------------------------------
# ROIs2follow
rois2follow<-read.table(paste0(input_dir,"all/reading/destrieux/summary/tables/ROImeasures2follow.txt")) %>%
  rename(model=V1)
#---------------------------------------------------
# define runs to compare
runs2compare=list.dirs(input_dir,recursive = FALSE)[grep("old|twins|comparison",list.dirs(input_dir,recursive = FALSE),invert=TRUE)] %>% gsub(input_dir,"",.) %>% gsub("\\/","",.)
# get paths
paths2compare=list.dirs(input_dir)
paths2compare=paths2compare[grep("old|twins|comparison",paths2compare,invert=TRUE)]
paths2compare=paths2compare[grep(paste(pheno,"/",atlas,"*.*","tables",sep=""),paths2compare)]
paths2compare=paths2compare[grep("summary",paths2compare)]
#---------------------------------------------------
for (r in runs2compare){
  ps<-paths2compare[grep(paste0("/",r,"/"),paths2compare)]
  
  if (length(ps)>0){
  for (p in ps){
    ts<-list.files(p,pattern="_t_table.csv")
    ts<-ts[grep("sig|rois2followup",ts,invert=TRUE)]
    for (t in ts){
      tmp<-read.csv(paste(p,t,sep="/"))
      tmp$subset<-r
      tmp$table<-t
      if (exists("d_all")){
        d_all<-merge(d_all,tmp,all=TRUE)
      } else {
        d_all<-tmp
      }
      rm(tmp)
    }
    rm(t,ts)
  }
    rm(p)
  }
  rm(ps)
}
rm(r)

#
d_all<-d_all %>%
  mutate(adj=factor(adj,levels=c("none","TotalAreaLH","MeanThicknessLH"))) %>% filter(!is.na(adj))
# filter out meanCT and totalCSA measures when adjusted (i.e. adjusted for themselves)
d_all<-d_all %>% filter(
  !(adj=="TotalAreaLH"&region=="total") &
  !(adj=="MeanThicknessLH"&region=="mean")
)
# rename labels to match text
levels(d_all$adj)<-c("Model 1","Model 2a\n(adj. for total CSA)","Model 2b\n(adj. for mean CT)")

# create additional columns for organization/visualization
d_all$model<-as.character(d_all$model)
d_all$model[is.na(d_all$model)]<-"baseline"
d_all$measure<-d_all$region<-NA
w<-which(d_all$model!="baseline")
if(length(w)>0){
  d_all$measure[w]<-strsplit(as.character(d_all$model[w]),"_") %>% sapply("[[",2)
  d_all$region[w]<-strsplit(as.character(d_all$model[w]),atlas) %>% sapply("[[",2) %>% gsub("^_","",.) %>% gsub("\\.rh|\\.lh","",.)
  
}
rm(w)

d_all<- d_all %>% mutate(roi2follow=if_else(model %in% rois2follow$model,region,""))
d_all$roi2follow[which(d_all$roi2follow=="")]<-NA

# from long to wide format, using subset for comparison
# long to wide
d_wide<- d_all %>% dplyr::select(Estimate,StdError,df,t.value,P,roi2follow,subset,model,measure,region,adj) %>% 
  pivot_wider(names_from=subset,values_from=c(Estimate,StdError,df,t.value,P))

d_wide<- d_wide %>% 
  mutate(label2show=if_else(is.na(roi2follow),1,2) %>% as.factor(),
         adjGlobal=if_else(adj!="Model 1",1,0) %>% as.factor())

# get correlations across ROIs, per model and per measure
d_wide_cors <- d_wide %>% group_by(measure,adj) %>%
  summarise(cor_Est_allVSallunrel=cor.test(Estimate_all,Estimate_all_unrelated)[["estimate"]] %>% round(.,digits=3),
            p_Est_allVSallunrel=cor.test(Estimate_all,Estimate_all_unrelated)[["p.value"]],
            cor_Est_allVSEURunrel=cor.test(Estimate_all,Estimate_eur6sd_unrelated)[["estimate"]] %>% round(.,digits=3),
            p_Est_allVSEURunrel=cor.test(Estimate_all,Estimate_eur6sd_unrelated)[["p.value"]],
            # t-values
            cor_T_allVSallunrel=cor.test(t.value_all,t.value_all_unrelated)[["estimate"]] %>% round(.,digits=3),
            p_T_allVSallunrel=cor.test(t.value_all,t.value_all_unrelated)[["p.value"]],
            cor_T_allVSEURunrel=cor.test(t.value_all,t.value_eur6sd_unrelated)[["estimate"]] %>% round(.,digits=3),
            p_T_allVSEURunrel=cor.test(t.value_all,t.value_eur6sd_unrelated)[["p.value"]]  
            )
d_cors <- d_wide_cors %>% 
  pivot_longer(all_of(c(contains("cor_"),contains("p_"))),names_to="tmp",values_to="value") %>%
  mutate(col=strsplit(tmp,"_") %>% sapply("[[",1),
         stat=strsplit(tmp,"_") %>% sapply("[[",2) %>%
           gsub("^T$","t.value",.) %>%
           gsub("Est","Estimate",.) ,
         comparison=strsplit(tmp,"_") %>% sapply("[[",3),
         sub1=strsplit(comparison,"VS") %>% sapply("[[",1) %>%
           gsub("unrel"," unrelated",.) %>% 
           gsub("EUR", " Eur6sd",.) %>%
           gsub("^ ","",.) %>% gsub(" ","_",.),
         sub2=strsplit(comparison,"VS") %>% sapply("[[",2) %>%
           gsub("unrel"," unrelated",.) %>% 
           gsub("EUR", " eur6sd",.) %>%
           gsub("^ ","",.) %>% gsub(" ","_",.)
         ) %>%
  dplyr::select(-tmp) %>%
  pivot_wider(names_from=col,values_from=value) %>%
  mutate(corLab=paste0("r=",round(cor,digits=2),"; p=",format(p,scientific=TRUE,digits=2)))

#--------------
# visualize
#--------------
s1="all"
s2="all_unrelated"
s3="eur6sd_unrelated"
n1=9013
n2=7502
n3=4080


#------ beta estimates --------
val="Estimate"
val_lims<-range(d_all[,val])
val_lims<-c(val_lims[1]-1*mean(val_lims),val_lims[2]+1*mean(val_lims))
# s1 vs s2
estimates_comparison_plot1 <- d_wide %>%
  mutate(x=get(paste0(val,"_",s1)),
         y=get(paste0(val,"_",s2)),
         se.x=get(paste0("StdError","_",s1)),
         se.y=get(paste0("StdError","_",s2))
  ) %>%
  ggplot(aes(x=x,y=y,color=measure,label=roi2follow,alpha=label2show,shape=adjGlobal)) + 
  geom_abline(intercept = 0,alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.5) +
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.5) +
  geom_text_repel(min.segment.length = 0,
                  box.padding = 0.8,
                  point.size=5, max.overlaps = Inf,
                  # direction="y",
                  segment.color="grey80") +
  geom_errorbar(aes(ymin=y-se.y,ymax=y+se.y)) + 
  geom_errorbarh(aes(xmin=x-se.x,xmax=x+se.x)) + 
  geom_point(size=3) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s2) %>% filter(measure=="area")), 
            aes(x=Inf, y=val_lims[1]+0.5*mean(val_lims), 
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s2) %>% filter(measure=="thick")), 
            aes(x=Inf, y=val_lims[1],
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  facet_grid(~ adj ) +
  scale_color_manual(values=brewer.pal(12, "Paired")[c(4,8)],
                     labels=c("CSA","CT")) +
  scale_shape_manual(values=c(19,17),
                     labels=c("Not adjusted for global measures","Adjusted for global measures")) +
  labs(
    # title=paste0(gsub("_"," ",s1) %>% simpleCap, " vs ", gsub("_"," ",s2) %>% simpleCap),
    x= paste0(gsub("_"," ",s1) %>% simpleCap," (N=",n1,")"),
    y= paste0(gsub("_"," ",s2) %>% simpleCap," (N=",n2,")"),
    color="",shape="") +
  lims(x=val_lims,y=val_lims) +
  guides(alpha="none") +
  theme_minimal_grid(font_size = 16) +
  panel_border() +
  theme(legend.position="none") +
  NULL

# s1 vs s3
estimates_comparison_plot2 <- d_wide %>%
  mutate(x=get(paste0(val,"_",s1)),
         y=get(paste0(val,"_",s3)),
         se.x=get(paste0("StdError","_",s1)),
         se.y=get(paste0("StdError","_",s3))
  ) %>%
  ggplot(aes(x=x,y=y,color=measure,label=roi2follow,alpha=label2show,shape=adjGlobal)) + 
  geom_abline(intercept = 0,alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.5) +
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.5) +
  geom_text_repel(min.segment.length = 0,
                  box.padding = 0.8,
                  point.size=5, max.overlaps = Inf,
                  # direction="y",
                  segment.color="grey80") +
  geom_errorbar(aes(ymin=y-se.y,ymax=y+se.y)) + 
  geom_errorbarh(aes(xmin=x-se.x,xmax=x+se.x)) + 
  geom_point(size=3) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s3) %>% filter(measure=="area")), 
            aes(x=Inf, y=val_lims[1]+0.5*mean(val_lims), 
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s3) %>% filter(measure=="thick")), 
            aes(x=Inf, y=val_lims[1],
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  facet_grid(~ adj ) +
  scale_color_manual(values=brewer.pal(12, "Paired")[c(4,8)],
                     labels=c("CSA","CT")) +
  scale_shape_manual(values=c(19,17),
                     labels=c("Not adjusted for global measures","Adjusted for global measures")) +
  labs(
    # title=paste0(gsub("_"," ",s1) %>% simpleCap, " vs ", gsub("_"," ",s3) %>% simpleCap),
    x= paste0(gsub("_"," ",s1) %>% simpleCap," (N=",n1,")"),
    y= paste0(gsub("_"," ",s3) %>% simpleCap," (N=",n3,")"),
    color="",shape="") +
  lims(x=val_lims,y=val_lims) +
  guides(alpha="none") +
  theme_minimal_grid(font_size = 16) +
  panel_border() +
  theme(legend.position="none") +
  NULL
# combine both plots
title <- ggdraw() + draw_label(
  paste0(gsub("\\.","-",val) %>% simpleCap, " comparison across subsets"),
  fontface='bold',size=14)

comb<-plot_grid(
  estimates_comparison_plot1,
  estimates_comparison_plot2,
  ncol=1,
  labels=c("a","b") )
leg<-get_legend(estimates_comparison_plot1 + theme(legend.position="bottom"))

estimates_comparison_plot<-plot_grid(title,comb,
                                     leg, ncol=1,rel_heights = c(0.1,2,0.1)
                                     )
estimates_comparison_plot %>% ggsave(.,file=paste0("Comparison_",val,".png"),width = 10,height = 10)
rm(title,comb,leg)
#--------------
val="t.value"
val_lims<-range(d_all[,val])
val_lims<-c(val_lims[1]-1*mean(val_lims),val_lims[2]+1*mean(val_lims))
# s1 vs s2
tval_comparison_plot1 <- d_wide %>%
  mutate(x=get(paste0(val,"_",s1)),
         y=get(paste0(val,"_",s2)),
         se.x=get(paste0("StdError","_",s1)),
         se.y=get(paste0("StdError","_",s2)),
         label2show=if_else(is.na(roi2follow),1,2)
  ) %>%
  ggplot(aes(x=x,y=y,color=measure,label=roi2follow,alpha=label2show,shape=adjGlobal)) + 
  geom_abline(intercept = 0,alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.5) +
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.5) +
  geom_text_repel(min.segment.length = 0,
                  box.padding = 0.8,
                  point.size=5, max.overlaps = Inf,
                  # direction="y",
                  segment.color="grey80") +
  geom_point(size=3) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s2) %>% filter(measure=="area")), 
            aes(x=Inf, y=val_lims[1]+0.5*mean(val_lims), 
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s2) %>% filter(measure=="thick")), 
            aes(x=Inf, y=val_lims[1],
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  facet_grid(~ adj ) +
  scale_color_manual(values=brewer.pal(12, "Paired")[c(4,8)],
                     labels=c("CSA","CT")) +
  scale_shape_manual(values=c(19,17),
                     labels=c("Not adjusted for global measures","Adjusted for global measures")) +
  labs(
    # title=paste0(gsub("_"," ",s1) %>% simpleCap, " vs ", gsub("_"," ",s2) %>% simpleCap),
    x= paste0(gsub("_"," ",s1) %>% simpleCap," (N=",n1,")"),
    y= paste0(gsub("_"," ",s2) %>% simpleCap," (N=",n2,")"),
    color="",shape="") +
  lims(x=val_lims,y=val_lims) +
  guides(alpha="none") +
  theme_minimal_grid(font_size = 16) +
  panel_border() +
  theme(legend.position="none") +
  NULL

# s1 vs s3
tval_comparison_plot2 <- d_wide %>%
  mutate(x=get(paste0(val,"_",s1)),
         y=get(paste0(val,"_",s3)),
         se.x=get(paste0("StdError","_",s1)),
         se.y=get(paste0("StdError","_",s3))
  ) %>%
  ggplot(aes(x=x,y=y,color=measure,label=roi2follow,alpha=label2show,shape=adjGlobal)) + 
  geom_abline(intercept = 0,alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.5) +
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.5) +
  geom_text_repel(min.segment.length = 0,
                  box.padding = 0.8,
                  point.size=5, max.overlaps = Inf,
                  # direction="y",
                  segment.color="grey80") +
  geom_point(size=3) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s3) %>% filter(measure=="area")), 
            aes(x=Inf, y=val_lims[1]+0.5*mean(val_lims), 
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  geom_text(data=(d_cors %>% filter(stat==val&sub1==s1&sub2==s3) %>% filter(measure=="thick")), 
            aes(x=Inf, y=val_lims[1],
                label=corLab, color=measure),
            hjust=1,
            inherit.aes = FALSE) +
  facet_grid(~ adj ) +
  scale_color_manual(values=brewer.pal(12, "Paired")[c(4,8)],
                     labels=c("CSA","CT")) +
  scale_shape_manual(values=c(19,17),
                     labels=c("Not adjusted for global measures","Adjusted for global measures")) +
  labs(
    # title=paste0(gsub("_"," ",s1) %>% simpleCap, " vs ", gsub("_"," ",s3) %>% simpleCap),
    x= paste0(gsub("_"," ",s1) %>% simpleCap," (N=",n1,")"),
    y= paste0(gsub("_"," ",s3) %>% simpleCap," (N=",n3,")"),
    color="",shape="") +
  lims(x=val_lims,y=val_lims) +
  guides(alpha="none") +
  theme_minimal_grid(font_size = 16) +
  panel_border() +
  theme(legend.position="none") +
  NULL
# combine both plots
title <- ggdraw() + draw_label(
  paste0(gsub("\\.","-",val) %>% simpleCap, " comparison across subsets"),
  fontface='bold',size=14)

comb<-plot_grid(
  tval_comparison_plot1,
  tval_comparison_plot2,
  ncol=1,
  labels=c("a","b") )
leg<-get_legend(tval_comparison_plot1 + theme(legend.position="bottom"))

tval_comparison_plot<-plot_grid(title,comb,
                                     leg, ncol=1,rel_heights = c(0.1,2,0.1)
)
tval_comparison_plot %>% ggsave(.,file=paste0("Comparison_",val,".png"),width = 10,height =10)
rm(title,comb,leg)

#--------------
