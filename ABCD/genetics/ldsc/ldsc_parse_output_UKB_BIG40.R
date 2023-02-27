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
primary_dir=paste0(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/")
genlang_dir=paste0(dir,"/coeduca/data/working_data/GenLang/projects/GWASMA_freeze_v1/data/")
howe_dir=paste0(dir,"resources/datasets/GWAS_sumstats/downloaded_data/ieu-openGWAS/")

ldsc_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/genotyping/ldsc/")
# ldsc_dir=paste0(dir,"resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/output/")

# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste0(dir,"general_scripts/helper_functions.R"))
source(paste0(dir,"general_scripts/genotyping/ldsc/helper_functions_ldsc.R"))
# read ggseg images and color table, for consistency
pheno_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/")
load(paste0(pheno_dir,"ggseg_brain_ROIs.Rdata"))
braincolor_codes<- braincolor_codes %>%
  mutate(region_name=Region  %>% gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% toupper() ) %>%
  select(label,color,measure,region_name,Region,hemi)

#----------------------------------------------------------------------
# read phenotype information
idps<-read.csv(paste0(primary_dir,"IDPs_summary.csv"),header=TRUE)
idps$pheno<-numeric_nchar(idps$Pheno)
idps$region<- idps$region %>% gsub("SURFACE","",.)
idps<-merge(idps,braincolor_codes,by.x=c("region","measure"),by.y=c("region_name","measure"),all.x=TRUE)
# lateralized measures
idps_lat<-read.table(paste0(primary_dir,"IDPs_lat_summary.txt"),header=TRUE) %>% filter(is.na(NotLat)) %>% select(-NotLat)
idps_lat$L[!is.na(idps_lat$L)]<-numeric_nchar(idps_lat$L[!is.na(idps_lat$L)])
idps_lat$R[!is.na(idps_lat$R)]<-numeric_nchar(idps_lat$R[!is.na(idps_lat$R)])
idps_lat<-merge(idps_lat,idps[,c("pheno","region","lobe","label","color","measure")],by.x=c("L","measure"),by.y=c("pheno","measure"))

# genlang runs info
runs_genlang<-read.csv(paste0(genlang_dir,"GenLang_freeze1_runs_overview.csv"))

#----------------------------------------------------------------------
# read summary result form LDSC run and combine with IDP info
#----------------------------------------------------------------------
## heritability
h2<-read_ldsc_h2(h2file=paste0(ldsc_dir,"summary_h2_ldsc.table"))
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
rg<-read_ldsc_rg(rgfile=paste0(ldsc_dir,"summary_LRrg_ldsc.table"))
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
rg_atlas<-read_ldsc_rg(rgfile=paste0(ldsc_dir,"summary_rg_ldsc.table"))
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
rg_cogn<-read_ldsc_rg(rgfile=paste0(ldsc_dir,"summary_cognitive_rg_ldsc.table"))
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
  mutate(pheno= (paste0(name,"_",subset,"_","STERR_GCOFF") %>% gsub("_all_","_combined_",.) ) ) %>% 
  select(pheno,name,subset,phenotype,measure)

rg_read<-read_ldsc_rg(rgfile=paste0(ldsc_dir,"summary_readingGenLang_rg_ldsc.table")) 
rg_read$p2<-rg_read$p2 %>% strsplit(.,"GWASMA_freeze_dald") %>% sapply("[[",2)
rg_read<-merge(rg_read,runs_genlang,by.x=c("p2"),by.y=c("pheno"),all.x=TRUE) # %>% filter(measure=="RT") 
rg_read<-merge(rg_read,idps,by.x=c("p1"),by.y=c("pheno"),all.x=TRUE,suffixes=c(".GenLang","")) # %>% filter(hemisphere=="L"|is.na(hemisphere))

#----------------------------------------------------------------------
## rg with sibGWAS from Howe et al. (2022)
## traits: Education Years and Cognitive Function
howe2022_traits<-read.table(paste0(howe_dir,"Howe2022_GWASlist2download.txt"),header=TRUE) %>%
  mutate(study=paste0(strsplit(author," ") %>% sapply("[[",1),year),
         measure.Howe2022=paste(trait,note)) %>%
  dplyr::select(study,id,trait,note,measure.Howe2022,sample_size)
rg_sibs<-read_ldsc_rg(rgfile=paste0(ldsc_dir,"summary_Howe2022_rg_ldsc.table"))
rg_sibs<-rg_sibs %>%
  mutate(p2= p2 %>% strsplit(.,"ieu-b-") %>% sapply("[[",2) %>% gsub(".tsv.gz","",.),
         id=paste0("ieu-b-",p2)) %>%
  merge(.,howe2022_traits,by=c("id"))
rg_sibs<-merge(rg_sibs,idps,by.x=c("p1"),by.y=c("pheno"),all.x=TRUE,suffixes=c(".Howe2022","")) # %>% filter(hemisphere=="L"|is.na(hemisphere))

#----------------------------------------------------------------------
write.csv(h2,paste0(ldsc_dir,"ldsc_h2_UKB_BIG40.csv"),row.names=FALSE)
write.csv(rg,paste0(ldsc_dir,"ldsc_rgLR_UKB_BIG40.csv"),row.names=FALSE)
write.csv(rg_atlas,paste0(ldsc_dir,"ldsc_rg_ROInGlobal_UKB_BIG40.csv"),row.names=FALSE)
write.csv(rg_cogn,paste0(ldsc_dir,"ldsc_rg_cognitive_UKB_BIG40.csv"),row.names=FALSE)
write.csv(rg_read,paste0(ldsc_dir,"ldsc_rg_readingGenLang_UKB_BIG40.csv"),row.names=FALSE)
write.csv(rg_sibs,paste0(ldsc_dir,"ldsc_rg_sibGWASHowe2022_UKB_BIG40.csv"),row.names=FALSE)

#--------------------------------------------