library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggseg)
library(ggsegDesterieux)

colors_brain<-function(data,atlas="destrieux",p_col="p1"){
  data$p_col=data[,p_col]
  data<- data %>% select(region,measure,hemisphere,p_col) %>% 
    mutate(region=tolower(region) %>% gsub(".lh|.rh","",.),measure=toupper(measure)) %>%
    distinct()
  # define colors - all
  colors<-c(brewer.pal(length(unique(data$region[which(data$measure=="THICKNESS")])), "Set1"),
            brewer.pal(length(unique(data$region[which(data$measure=="AREA")]))+1, "Set2")[-2])
  
  # manually replace two colors, to avoid bright yellow...
  colors[which(colors=="#FFFF33")]<-"#B3B3B3"
  
  # associate color to region
  if(atlas=="destrieux"){
    target_labels<-data.frame(target=ggseg(atlas=desterieux)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
  } else if (atlas=="desikan"){
    target_labels<-data.frame(target=ggseg(atlas=dkt)$data$label %>% unique(),stringsAsFactors = FALSE)  %>% mutate (matchLabel=tolower(target) %>% gsub("\\.|-","_",.))
  }
  
  rois_summary <- data %>% 
    mutate(hemi=if_else(str_detect(p_col,"lh"),"lh",if_else(str_detect(p_col,"rh"),"rh",""))) %>%
    filter(hemi!="") %>%
    select(region,measure,hemi,hemisphere) %>%
    mutate(matchLabel=paste(hemi,region,sep="_") %>% 
             gsub("area_|thick_|smri_|cort.","",.) %>% 
             gsub(paste0(atlas,"_"),"",.) %>%
             gsub("\\.|-","_",.) ) %>%
    mutate(hemi=if_else(hemi=="lh","left",if_else(hemi=="rh","right",""))) %>% 
    distinct() %>% arrange(measure, matchLabel,hemi)
  colorcode<-data.frame(region=unique(rois_summary$region),color=colors)
  rois_summary<-merge(rois_summary,colorcode,by="region")
  
  rois_summary <- merge(rois_summary,target_labels,by="matchLabel",stringsAsFactors=TRUE,all.x=TRUE) %>% 
    mutate(label=as.factor(target)) %>% select(label,color,measure,region,hemi,hemisphere) %>% 
    rename(Region=region)  %>% arrange(measure,Region)
  
  # define order based on labels, which will be constant regardless of original region names
  rois_summary$label2<-as.character(rois_summary$label) %>% gsub("rh_","lh_",.) # because same color to each reg. regardless of hemi
  rois_summary$label2[is.na(rois_summary$label2)]<-"NA"
  rois_summary$label2<-factor(rois_summary$label2,
                             # specify order of levels for consistency of colors, and order in plots
                             levels=c("NA","lh_G_parietal_sup","lh_G_temp_sup-Lateral",
                                      "lh_G_front_middle","lh_G_front_sup","lh_G_postcentral",
                                      "lh_G_oc-temp_lat-fusifor","lh_G_oc-temp_med-Lingual","lh_G_and_S_occipital_inf"))
  rois_summary<-rois_summary %>% arrange(label2)
  rois_summary$Region<-factor(rois_summary$Region,
                              # specify order of levels for consistency of colors, and order in plots
                              levels=c(unique(rois_summary$Region))) 
  rois_summary<-rois_summary%>%arrange(Region)
  return(rois_summary)
}

