library(dplyr)

#
args <- commandArgs(TRUE)
# args<-c(
#   # info file
#   # "F:/projects/resources/datasets/ABCD/data/working_data/phenotypes/ABCD_release2.0.1_sex_ancestry.txt",
#   "F:/projects/resources/datasets/ABCD/data/primary_data/ABCDStudyNDAv3/acspsw03.txt",
#   # fam file
#   "F:/projects/resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/ABCD_release_3.0_QCed.updated_genders.fam"
#   )


# args<-c(
#   # info file
#   # "F:/projects/resources/datasets/ABCD/data/working_data/phenotypes/ABCD_release2.0.1_sex_ancestry.txt",
#   "F:/projects/resources/datasets/ABCD/data/primary_data/ABCDStudyNDAv3/acspsw03.txt",
#   # fam file
#   "F:/projects/resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/backup/intermediate_files/ABCD_release_3.0_QCed.updated_genders.fam"
#   )

# extract genotyping dir, just in case from file2
geno_dir<-sapply(strsplit(args[2],"/ABCD_"),"[[",1) %>% unique()


#---------------------------------------#
# READ DATA                             #
#---------------------------------------#
info<-read.table(args[1],header=TRUE,stringsAsFactors = FALSE) %>% 
                filter(collection_id!="collection_id") %>% 
                filter(eventname=="baseline_year_1_arm_1") %>%
                rename(IID=src_subject_id) %>%
                mutate(AFR=genetic_af_african %>% as.numeric(),
                       EUR=genetic_af_european %>% as.numeric(),
                       EAS=genetic_af_east_asian %>% as.numeric(),
                       AMR=genetic_af_american %>% as.numeric())
info_description<-read.table(args[1],header=TRUE,stringsAsFactors = FALSE,nrow=1)
fam<-read.table(args[2],header=FALSE)
colnames(fam)<-c("FID","IID","pID","mID","sex","pheno")

# combine
info_f<-merge(info,fam,by=c("IID"),suffixes=c(".info",".fam"),all.y=TRUE)

#---------------------------------------#
# Define ancestry groups                 #
#---------------------------------------#
# check for each level of ethnicity, what is the proportion of each of the four ancestry groups
summary_group<-lapply(unique(info_f$race_ethnicity), function(x){
  info_f %>% filter(race_ethnicity==x) %>% select("EUR","AFR","EAS","AMR") %>% summary() 
})

apply(info[,c("EUR","AFR","EAS","AMR")],2,summary) %>% as.data.frame()
apply(info[,c("EUR","AFR","EAS","AMR")],2,function(x) {table(is.na(x))})
# select individuals with a 90% of one ancestry (or unkown), to be analyzed as separate groups
# include also if genetic ancestry is not defined
thr=0.9
lapply(c("EUR","AFR","EAS","AMR"),function(x){
  w1<-which(info_f[,x]>thr)
  d<-info_f[w1,c("FID","IID",x)]
  write.table(d %>% select(FID,IID),file=paste("ancestry_",thr,x,".samples",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
  return(dim(d))
})

#---------------------------------------#
