# clean workspace
rm(list=ls())
# libraries and custom functions
library(dplyr)
library(tidyr)

#---------------------------------------------------
options(stringsAsFactors = FALSE)
root="UKB_BIG40"
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")
genlang_dir=paste(dir,"/coeduca/data/working_data/GenLang/projects/GWASMA_freeze_v1/data/",sep="")
# ldsc_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/ldsc/",sep="")
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/",sep="")
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
source(paste(dir,"general_scripts/genotyping/ldsc/helper_functions_ldsc.R",sep=""))
# read ggseg images and color table, for consistency
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/",sep="")
load(paste(pheno_dir,"ggseg_brain_ROIs.Rdata",sep=""))
braincolor_codes<- braincolor_codes %>%
  mutate(region_name=Region  %>% gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% toupper() ) %>%
  select(label,color,measure,region_name,Region,hemi)

#----------------------------------------------------------------------
# read phenotype information
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE)
idps$pheno<-numeric_nchar(idps$Pheno)
idps$region<- idps$region %>% gsub("SURFACE","",.)
idps<-merge(idps,braincolor_codes,by.x=c("region","measure"),by.y=c("region_name","measure"),all.x=TRUE)
# lateralized measures
idps_lat<-read.table(paste(primary_dir,"IDPs_lat_summary.txt",sep=""),header=TRUE) %>% filter(is.na(NotLat)) %>% select(-NotLat)
idps_lat$L[!is.na(idps_lat$L)]<-numeric_nchar(idps_lat$L[!is.na(idps_lat$L)])
idps_lat$R[!is.na(idps_lat$R)]<-numeric_nchar(idps_lat$R[!is.na(idps_lat$R)])
idps_lat<-merge(idps_lat,idps[,c("pheno","region","lobe","label","color","measure")],by.x=c("L","measure"),by.y=c("pheno","measure"))

# genlang runs info
runs_genlang<-read.csv(paste0(genlang_dir,"GenLang_freeze1_runs_overview.csv"))

#----------------------------------------------------------------------
# read summary result form LDSC run and combine with IDP info
#----------------------------------------------------------------------
## heritability
h2<-read_ldsc_h2(h2file=paste(ldsc_dir,"summary_h2_ldsc.table",sep=""))
# combine with IDP info and select columns to keep
h2<-merge(h2,idps,by=c("pheno"),all.x=TRUE)
h2<-h2 %>% select(IDP_short_name,pheno,parcellation,measure,color,lobe,region,hemisphere,h2,h2.SE,Intercept,Intercept.SE,Lambda.GC,h2.CI95upper,h2.CI95lower,N) %>%
  arrange(lobe,region,measure,parcellation)
h2$p<-pchisq((h2$h2/h2$h2.SE)^2,1,lower.tail = FALSE) %>% format(scientific=TRUE,digits=2) %>% as.numeric()
# h2$h2<-h2$h2 %>% format(digits=3)  %>% as.numeric()
# h2$h2.SE<-h2$h2.SE %>% format(digits=3)  %>% as.numeric()
h2$measure<-as.factor(h2$measure)

#----------------------------------------------------------------------
## genetic correlation between left and right (for lateral measures)
## left and right
rg<-read_ldsc_rg(rgfile=paste(ldsc_dir,"summary_LRrg_ldsc.table",sep=""))
rg<-rg%>% rename(R=p1,L=p2)
# combine with IDP info and select columns to keep
rg<-merge(idps_lat,rg,by.x=c("L","R","IDP_short_name_lat"),by.y=c("L","R","file"),all.y=TRUE)
rg<-rg %>% select(IDP_short_name_lat,parcellation,measure,color,lobe,region,L,R,rg,rg.SE,rg.CI95upper,rg.CI95lower,p,gcov_int,gcov_int_se) %>%
  arrange(lobe,region,measure,parcellation) %>% filter(parcellation!="wg")
rg$p2<-pchisq(((1-rg$rg)/rg$rg.SE)^2,1,lower.tail = FALSE) %>% format(scientific=TRUE,digits=2) %>% as.numeric() # statistical significance for different to 1
# rg$rg<-rg$rg %>% format(digits=3)  %>% as.numeric()
# rg$rg.SE<-rg$rg.SE %>% format(digits=2)  %>% as.numeric()
rg$measure<-as.factor(rg$measure)

#----------------------------------------------------------------------
# ## rg within parcellation and measure type (left and right)
rois2follow <- h2$region %>% unique()
rg_atlas<-read_ldsc_rg(rgfile=paste(ldsc_dir,"summary_rg_ldsc.table",sep=""))
# rg_atlas<-subset(rg_atlas,p1=="0648"|p1=="1020")
# combine with IDP info and select columns to keep
rg_atlas<-merge(rg_atlas,idps,by.x=c("p1"),by.y=c("pheno"),all.x=TRUE) %>%
  mutate(region=gsub("GLOBALMEANTHICKNESS","MEAN",region))
rg_atlas<-merge(rg_atlas,idps,by.x=c("p2"),
                by.y=c("pheno"),
                suffixes=c(".p1",".p2"),all.x=TRUE) %>% filter(!is.na(rg))
rg_atlas<-rg_atlas %>% filter(region.p2 %in% rois2follow)
# only between same measure types
rg_atlas<-rg_atlas[which(rg_atlas$measure.p1==rg_atlas$measure.p2),]
rg_atlas<-rg_atlas %>% rename(measure=measure.p1) %>% select(-measure.p2)

rg_atlas<-rg_atlas %>% select(IDP_short_name.p1,IDP_short_name.p2,measure,region.p1,region.p2,p1,p2,
                              rg,rg.SE,rg.CI95upper,rg.CI95lower,p,gcov_int,gcov_int_se)  %>%
  arrange(measure)
rg_atlas$rg<-rg_atlas$rg %>% format(digits=1)  %>% as.numeric()
rg_atlas$rg.SE<-rg_atlas$rg.SE %>% format(digits=2)  %>% as.numeric()
rg_atlas$measure<-as.factor(rg_atlas$measure)

#----------------------------------------------------------------------
## rg with cognitive measures
rg_cogn<-read_ldsc_rg(rgfile=paste(ldsc_dir,"summary_cognitive_rg_ldsc.table",sep=""))
rg_cogn<-merge(rg_cogn,idps,by.x=c("p1"),by.y=c("pheno"),all.x=TRUE)  %>% 
  filter(p2!="Savage2018_IQ")
rg_cogn<-rg_cogn %>% select(IDP_short_name,parcellation,hemisphere,measure,color,region,
                            p1,p2,rg,rg.SE,rg.CI95upper,rg.CI95lower,p,gcov_int,gcov_int_se)  %>%
  arrange(measure,parcellation)
rg_cogn$rg<-rg_cogn$rg %>% format(digits=1)  %>% as.numeric()
rg_cogn$rg.SE<-rg_cogn$rg.SE %>% format(digits=2)  %>% as.numeric()
rg_cogn$measure<-as.factor(rg_cogn$measure)
rg_cogn$cognitive<-NA
rg_cogn$cognitive[rg_cogn$p2=="Lee2018_EA3"]<-"Educational Attainement"
rg_cogn$cognitive[rg_cogn$p2=="Lee2018_CP"]<-"Cognitive Performance"
rg_cogn$cognitive[rg_cogn$p2=="DDyslexia_Gialluisi2020"]<-"Developmental Dyslexia"

#----------------------------------------------------------------------
## rg with reading measures
runs_genlang <- runs_genlang %>% 
  mutate(pheno= (paste(name,"_",subset,"_","STERR_GCOFF",sep="") %>% gsub("_all_","_combined_",.) ) ) %>% 
  select(pheno,name,subset,phenotype,measure)

rg_read<-read_ldsc_rg(rgfile=paste(ldsc_dir,"summary_readingGenLang_rg_ldsc.table",sep="")) 
rg_read$p2<-rg_read$p2 %>% strsplit(.,"GWASMA_freeze_dald") %>% sapply("[[",2)
rg_read<-merge(rg_read,runs_genlang,by.x=c("p2"),by.y=c("pheno"),all.x=TRUE) %>% filter(measure=="RT") 
rg_read<-merge(rg_read,idps,by.x=c("p1"),by.y=c("pheno"),all.x=TRUE) # %>% filter(hemisphere=="L"|is.na(hemisphere))

#----------------------------------------------------------------------
write.csv(h2,paste(ldsc_dir,"ldsc_h2_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rg,paste(ldsc_dir,"ldsc_rgLR_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rg_atlas,paste(ldsc_dir,"ldsc_rg_ROInGlobal_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rg_cogn,paste(ldsc_dir,"ldsc_rg_cognitive_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rg_read,paste(ldsc_dir,"ldsc_rg_readingGenLang_UKB_BIG40.csv",sep=""),row.names=FALSE)

#--------------------------------------------