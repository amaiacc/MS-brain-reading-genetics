# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/reading.config" # parameters for this run
hemi="lh"

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
library("ggseg")
library("ggsegExtra")
# ggseg::install_atlases()
library(ggsegDesterieux)
library(ggsegDefaultExtra)
#
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())

# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
# source config file to define parameters for this run
source(config_file)
#---------------------
if(atlas=="destrieux"){
  target_labels<-data.frame(target=ggseg(atlas=desterieux)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
} else if (atlas=="desikan"){
  target_labels<-data.frame(target=ggseg(atlas=dkt)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
}

base_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",atlas,"/",sep="")
setwd(base_dir)
#---------------------

# get ROIs2follow, for main figures
rois2follow<-read.table("summary/tables/ROImeasures2follow.txt") %>% rename(model=V1)

# baseline models
ttables<-c( paste(paste("baseline","/tables/",sep=""),
                  list.files(paste("baseline","/tables/",sep=""),pattern="t_table"),sep="/"))

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
  t<-merge(t,target_labels,by="matchLabel",all.x=TRUE) %>% mutate(label=as.factor(target)) %>% select(-matchLabel,-target) 
  t<-t %>% rename(P= Pr...t..,StdError=Std..Error) %>% mutate(minuslog10P=(-log10(P)))
  t$label<-as.factor(t$label)
  # t<- t %>% filter(region!="total" & region!="mean")
  t$Pfdr<-p.adjust(t$P,method = "fdr")
  t$Pbonf<-p.adjust(t$P,method = "bonferroni")
  t$adj<-c
  # combine data
  if(exists("d")){d<-rbind(d,t)} else {d<-t}
  
  # make figures
  for (m in c(unique(t$measure))){
    
    titlep<-paste(dv,"~\n",c,"\n + ",m, " ROI",sep="")  # "(",hemi,")"
    t2<-t %>% filter(measure==m &hemi==hemi&region!="total"&region!="mean")
    t2$rois2follow<-FALSE
    t2$rois2follow[t2$model %in% rois2follow$model ]<-TRUE

    p<- t2 %>% select(label,t.value,P,model) %>% brain_join(desterieux) %>% filter(hemi=="left") %>%
      ggplot() + 
      geom_sf(aes(fill=t.value),colour="black") +
      scale_fill_gradient2(na.value="transparent",name="T-value",limit=c(-10,10)) +
      theme_void() +
      labs(caption=titlep) +
      NULL
    
    pfdr<- t2 %>% filter(Pfdr<0.05) %>% select(label,t.value,P,model) %>% brain_join(desterieux) %>% filter(hemi=="left") %>%
      ggplot() + 
      geom_sf(aes(fill=t.value),colour="black") +
      scale_fill_gradient2(na.value="transparent",name="T-value",limit=c(-10,10)) +
      theme_void() +
      labs(caption=titlep) +
      NULL
    
    pROIs2follow <- t2 %>% filter(rois2follow==TRUE) %>% select(label,t.value,P,model) %>% brain_join(desterieux) %>% filter(hemi=="left") %>%
      ggplot() + 
      geom_sf(aes(fill=t.value),colour="black",size=1) +
      geom_sf_label(aes(label = ifelse(model %in% rois2follow$model, region, NA)),
                    alpha = .8,
                    show.legend = FALSE)  + 
      scale_fill_gradient2(na.value="transparent",name="T-value",limit=c(-10,10)) +
      theme_void() +
      labs(caption=titlep) +
      NULL
     
    assign(paste("plot_",m,"_m",i,sep=""),p)
    assign(paste("plot_",m,"_m",i,"_FDRsig",sep=""),pfdr)
    assign(paste("plot_",m,"_m",i,"_ROIs",sep=""), pROIs2follow)
    rm(p,pfdr,pROIs2follow)
    rm(t2)
  }
  rm(m,titlep)
}
rm(t)


#---------------------
# save combined plots
working_dir=paste(base_dir,"baseline","/",sep="")
out_dir=paste(working_dir,"/figures/",sep="")
if(!dir.exists(out_dir)){dir.create(out_dir)}

# make figures combining different models
thick_list<-ls(pattern="plot_thick_m")
area_list<-ls(pattern="plot_area_m")
# extract a legend that is laid out horizontally
legend_t <- get_legend(  get(thick_list[1]) +
                           guides(color = guide_legend(nrow = 1)) +  theme(legend.position = "bottom") )
legend_a <- get_legend(  get(area_list[1]) +
                           guides(color = guide_legend(nrow = 1)) +  theme(legend.position = "bottom") )
#
n_t<-ceiling(length(thick_list[grep("sig",thick_list,invert=TRUE)])/4)
n_a<-ceiling(length(area_list[grep("sig",thick_list,invert=TRUE)])/4)


#---------------------
# main figures - focusing on results of interest
## thickness

thick_ROIs<-plot_grid( plot_grid(plot_thick_m0_ROIs  + theme(legend.position="none"),
                                   plot_thick_m0d_ROIs + theme(legend.position="none"),
                                   nrow=2,
                                   labels=c("a","b")),  legend_t , ncol = 1, rel_heights = c(1, .1))

# area
area_ROIs<-plot_grid( plot_grid(plot_area_m0_ROIs  + theme(legend.position="none"),
                                   plot_area_m0c_ROIs + theme(legend.position="none"),
                                   nrow=2,
                                  labels=c("a","b")),  legend_t , ncol = 1, rel_heights = c(1, .1))


ggsave(thick_ROIs,file=paste(out_dir,pheno,"_",atlas,"_","thick_ROIs_main",".pdf",sep=""))
ggsave(area_ROIs,file=paste(out_dir,pheno,"_",atlas,"_","area_ROIs_main",".pdf",sep=""))

#---------------------
# all figures

plot_thick_all<-plot_grid( plot_grid(plotlist =lapply(thick_list[grep("sig|ROI",thick_list,invert=TRUE)], 
                                                      function(x) {get(x) + theme(legend.position="none")}),
                                     nrow=n_t),  legend_t , ncol = 1, rel_heights = c(1, .1)) 

plot_thick_FDRsig<-plot_grid( plot_grid(plotlist =lapply(thick_list[grep("FDRsig",thick_list,invert=FALSE)], 
                                                         function(x) {get(x) + theme(legend.position="none")}),
                                        nrow=n_t),  legend_t , ncol = 1, rel_heights = c(1, .1))

plot_thick_ROIs<-plot_grid( plot_grid(plotlist =lapply(thick_list[grep("ROIs",thick_list,invert=FALSE)], 
                                                       function(x) {get(x) + theme(legend.position="none")}),
                                      nrow=n_t),  legend_t , ncol = 1, rel_heights = c(1, .1)) 

plot_area_all<-plot_grid( plot_grid(plotlist =lapply(area_list[grep("sig|ROI",area_list,invert=TRUE)], 
                                                     function(x) {get(x) + theme(legend.position="none")}),
                                    nrow=n_a),  legend_a , ncol = 1, rel_heights = c(1, .1)) 

plot_area_FDRsig<-plot_grid( plot_grid(plotlist =lapply(area_list[grep("FDRsig",area_list)], 
                                                        function(x) {get(x) + theme(legend.position="none")}),
                                       nrow=n_a),  legend_a , ncol = 1, rel_heights = c(1, .1)) 
plot_area_ROIs<-plot_grid( plot_grid(plotlist =lapply(area_list[grep("ROIs",area_list,invert=FALSE)], 
                                                       function(x) {get(x) + theme(legend.position="none")}),
                                      nrow=n_a),  legend_a, ncol = 1, rel_heights = c(1, .1))
# save combined plots
for (m in c("area","thick")){
  # pdf(file=paste(out_dir,pheno,"_",atlas,"_",m,"_main_FDRsig",".pdf",sep=""),width=8,height=8)
  # paste(m,"_FDRsig",sep="") %>% get() %>% print()
  # dev.off() 
  
  
  pdf(file=paste(out_dir,pheno,"_",atlas,"_",m,".pdf",sep=""),width=15,height=7)
  paste("plot_",m,"_all",sep="") %>% get() %>% print()
  dev.off() 
  # pdf(file=paste(out_dir,pheno,"_",atlas,"_",m,"_BFsig",".pdf",sep=""),width=15,height=7)
  # paste("plot_",m,"_BFsig",sep="") %>% get() %>% print()
  # dev.off()  
  pdf(file=paste(out_dir,pheno,"_",atlas,"_",m,"_FDRsig",".pdf",sep=""),width=15,height=7)
  paste("plot_",m,"_FDRsig",sep="") %>% get() %>% print()
  dev.off()  
  
  pdf(file=paste(out_dir,pheno,"_",atlas,"_",m,"_ROIs",".pdf",sep=""),width=15,height=7)
  paste("plot_",m,"_ROIs",sep="") %>% get() %>% print()
  dev.off()  
  
}
rm(m)

#---------------------
# clean intermediate
rm(list=thick_list)
rm(list=area_list)
rm(thick_list,area_list,ttables)
rm(n_t,legend_t,n_a,legend_a)
