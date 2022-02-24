# Define the following subsets of individuals, to be used in the ABCD analyses:
# - all
# - eur (6sd, 3sd)
# - unrelated
# - twins
#---------------------------------------------------

# clean workspace
rm(list=ls())
#---------------------------------------------------
# variables - currently hard coded
# genetic data release
root_name="ABCD_release_3.0_QCed" 
# define phenotypic data version
v= "ABCDv3_DEAP" # "ABCD_release2.0.1"

## uncomment to take variables as arguments
# args <- commandArgs(TRUE)
## genetic data release
# root_name<-args[1]
## define phenotypic data version
# v<-args[2]
#---------------------------------------------------

# define libraries
library(dplyr); library(tidyr)
#
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())

#---------------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/",sep="")
geno_clean_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/clean/",sep="")

working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",sep="")
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")

# get to working dir
setwd(working_dir)
#---------------------------------------------------
# read data
#---------------------------------------------------

#---------------------------------------------------
# Subset fam files, from genetic QC analysis: preImputation_batchQC_ABCD.sh
#---------------------------------------------------
all<-read.table(paste(geno_clean_dir,root_name,"_all.cleaned.fam",sep="")) %>% mutate(all=1)
all_unrelated<-read.table(paste(geno_clean_dir,root_name,"_all_unrelated.cleaned.fam",sep="")) %>%
  mutate(all_unrelated=1)
eur6sd<-read.table(paste(geno_clean_dir,root_name,"_EUR.EUR6sd.cleaned.fam",sep="")) %>%
  mutate(eur6sd=1)
eur6sd_unrelated<-read.table(paste(geno_clean_dir,root_name,"_EUR.EUR6sd_unrelated.cleaned.fam",sep="")) %>%
  mutate(eur6sd_unrelated=1)
eur3sd<-read.table(paste(geno_clean_dir,root_name,"_EUR.EUR3sd.cleaned.fam",sep="")) %>%
  mutate(eur3sd=1)
eur3sd_unrelated<-read.table(paste(geno_clean_dir,root_name,"_EUR.EUR3sd_unrelated.cleaned.fam",sep=""))%>%
  mutate(eur3sd_unrelated=1)
#---------------------------------------------------
## additional information, intermediate files
# cols selected in: subset_data_*.R
rel <-read.table(paste0(v,"_sex_ancestry.txt"),header=TRUE,sep="\t") %>%
  select(-AFR,-EUR,-EAS,-AMR)

table(rel$rel_relationship)
twins<-subset(rel,rel_relationship=="twin")

# relatedness, to define unrelated and twins sets from genetic data (based on pi-hat)
# pi_hats
pi_hat<-read.table(paste(geno_dir,"/relatedness/",root_name,"_pruned_ibd_pairs_pihat.txt",sep=""),header=TRUE) %>%
  filter(PI_HAT>0.12)

# related (randomly chosen from pairs 0.1875, and excluded for unrelated samples)
related<-read.table(paste(geno_dir,"/relatedness/",root_name,"_pruned_ibd_outliers.txt",sep=""),header=TRUE) %>%
  mutate(related=1)
related<-related[!duplicated(related),]
#---------------------------------------------------

#---------------------------------------------------
# combine files
subsets<-c("all","all_unrelated","eur6sd","eur6sd_unrelated","eur3sd","eur3sd_unrelated")
for (s in subsets){
  t<-get(s)
  if(!exists("subsets_geneticQC")){subsets_geneticQC<-t} else{ subsets_geneticQC<-merge(subsets_geneticQC,t,all=TRUE)}
  rm(t)
}
rm(s)
rm(all,all_unrelated,eur3sd,eur3sd_unrelated,eur6sd,eur6sd_unrelated)

colnames(subsets_geneticQC)<-c("FID","IID","pID","mID","sex","pheno",subsets)

subsets_geneticQC[,subsets][is.na(subsets_geneticQC[,subsets])]<-0

subsets_geneticQC<-subsets_geneticQC %>% 
                    select("FID","IID","sex",subsets) %>%
                    rename(subjectid=IID)

# add unrelated flag
subsets_geneticQC<-merge(subsets_geneticQC,related,by.x=c("subjectid","FID"),by.y=c("IID1","FID1"),all.x=TRUE)
subsets_geneticQC<- subsets_geneticQC %>% mutate(unrelated=if_else(is.na(related),1,0))

# add twin flag
subsets_geneticQC<-merge(subsets_geneticQC,twins[,c("subjectid","Zygosity")],by="subjectid",all.x=TRUE) %>%
                    mutate(twin=Zygosity %>% as.character()) 
table(subsets_geneticQC$twin,useNA="ifany")
# recode twin as 0/1
subsets_geneticQC$twin[!is.na(subsets_geneticQC$twin)]<-1
subsets_geneticQC$twin[is.na(subsets_geneticQC$twin)]<-0

# remove unnecessary columns
subsets_geneticQC<- subsets_geneticQC %>% select(-FID,-Zygosity,-related)

# save subset file
write.table(subsets_geneticQC,paste(working_dir,root_name,"_",v,"_subsets.table",sep=""),row.names = FALSE,quote=FALSE)
