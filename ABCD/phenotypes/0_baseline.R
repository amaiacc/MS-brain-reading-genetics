# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]

config_file="reading.config" # parameters for this run
# config_file="F:/projects//resources/datasets/ABCD/scripts/phenotypes/reading.config" # parameters for this run

#---------------------


#---------------------------------------------------
# define libraries
library(data.table)
library(dplyr); library(tidyr)
# library(gamm4) # 
library(moments) # skewness and kurtosis
#
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
#---------------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",sep="")
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")

# general functions to be used
source(paste(scripts_dir,"general_functions.R",sep=""))

# source config file to define parameters for this run
source(paste(scripts_dir,config_file,sep="")) # as parameter
# source(config_file) # if hard-coded

# get to working dir
setwd(working_dir)


#---------------------------------------------------
# get subset definitions, as defined in define_subsets.R
## and subset to selected individuals only (i.e. s)
id_subset<-read.table(paste0(gen_v,"_subsets.table"),header=TRUE) %>% select("subjectid",s)
id_subset<-id_subset[which(id_subset[,s]==1),]

# get PCs, from file defined in config
pcs<-fread(paste(geno_dir,pc_file,sep=""),header=FALSE) ## all samples
colnames(pcs)<-c("FID","subjectid",paste("PC",1:100,sep=""),"pheno")
pcs<-pcs %>% select(-pheno)

# read data
dat <- readRDS(paste0(v,"_baseline_sMRI.Rds")) 
# NROW(dat) # 11265
#---------------------------------------------------
# combine and subset data
d<-merge(dat,id_subset) # N=10455
d<-merge(d,pcs) # N=10454
d$rel_family_id<-as.factor(d$rel_family_id) 
#---------------------------------------------------

#---------------------------------------------------
# define variables and covariates
# defined in config file: 
# - covs_d ## covariates
# - cov11, cov12 and cov2 ## additional covariates
# - covs_int # interaction terms to consider
covs_d=c("sex","age","high.educ.bl","household.income.bl"
         # "married.bl", "race.4level", "hisp"
         )
covs_def<-c(covs_d,paste("PC",1:10,sep=""))


# independent variables to assess
# i.e., defines rois to be used in analysis
mri_info<-colnames(dat)[grep("mri_info*.*device",colnames(dat))]
total<-colnames(dat)[grep("total|mean",colnames(dat))] # measures are same for desikan
total<-total[grep(atlas,total)]
total<-total[grep("vol|thick|area",total)]
total_vol<-total[grep("vol*.*total$",total)]
total_vol23<-paste(total_vol,"2.3",sep="")
lh_area<-total[grep("area*.*lh",total)]
lh_thick<-total[grep("thick*.*.lh",total)]

## covariates to consider (for sensitivity)
covs_add<-c(total_vol,lh_area,lh_thick,cov11,cov12,cov2)

## all variables to consider
vars<-c(covs_def,dv,covs_add,
        total,
        mri_info,"rel_family_id","abcd_site") %>% unique() #ancestry_vars,

n_vars<-vars[grep("numeric",d[,vars] %>% sapply(.,class))]
c_vars<-vars[grep("character|factor",d[,vars] %>% sapply(.,class))]

# define rois
rois<-colnames(d)[grep(atlas,colnames(d))]
rois<-rois[grep("thick|area",rois)]
rois<-rois[grep("total",rois,invert=TRUE)]
rois_lh<-rois[grep("lh",rois)]
#---------------------------------------------------
# define complete data

## check missing data vars
apply(d[,vars],2,function(c) table(is.na(c)))

##
d0<-d[complete.cases(d[,c(vars,rois)]),] # N=9224

#---------------------------------------------------
# for each of the numeric variables, remove extreme outliers, i.e. if they are > thr*SD from the mean
n_vars2check<-n_vars[grep("PC",n_vars,invert=TRUE)]
d0[,n_vars] <- apply(d0[,n_vars],2,function(x){
  rm_outliers(x,thr=thr)
})

# ##if they are more extreme than the c(0.00025,0.99975) quantiles
# apply(d0[,n_vars],2,function(x,f=0.0005){
#   censor(x,fraction=f) %>% is.na() %>% sum()
# 
# })

# keep complete cases only
d0<-d0[complete.cases(d0[,vars]),] # N=9177

#---------------------------------------------------

#---------------------------------------------------
# ## QC
#---------------------------------------------------
# # double check that each family only has info from one scanner and abcd site
# counts_mri_per_family<-sapply(unique(d0$rel_family_id), function(x){
#   t<-d0 %>% filter(rel_family_id==x) %>% pull(mri_info_device.serial.number) %>% unique() %>% length()
#   return(t)
# })
# table(counts_mri_per_family)
# 
# counts_site_per_family<-sapply(unique(d0$rel_family_id), function(x){
#   t<-d0 %>% filter(rel_family_id==x) %>% pull(abcd_site) %>% unique() %>% length()
#   return(t)
# })
# table(counts_site_per_family)
# 
# fam_multiple_mris<-unique(d0$rel_family_id)[which(counts_mri_per_family>1)]
# fam_multiple_sites<-unique(d0$rel_family_id)[which(counts_site_per_family>1)]
# 
# lapply(unique(c(fam_multiple_mris,fam_multiple_sites)),function(f){
#   subset(d0,rel_family_id==f)[,c("subjectid","rel_family_id",mri_info,"abcd_site","PI_HAT")]
# })

#-----------------------------------------------------
## create output directories, if not there already
# each dir within path should be created first
if(!dir.exists(paste(working_dir,v,sep="/"))){dir.create(paste(working_dir,v,sep="/"))}
if(!dir.exists(paste(working_dir,v,s,sep="/"))){dir.create(paste(working_dir,v,s,sep="/"))}
if(!dir.exists(paste(working_dir,v,s,pheno,sep="/"))){dir.create(paste(working_dir,v,s,pheno,sep="/"))}
#
base_dir=paste(working_dir,v,s,pheno,paste("baseline",atlas,sep="_"),sep="/")
fig_dir=paste(base_dir,"/trim",trim_val,"/figures/",sep="")
if(!dir.exists(base_dir)){dir.create(base_dir)}
if(!dir.exists(paste(base_dir,"/trim",trim_val,"/",sep=""))){dir.create(paste(base_dir,"/trim",trim_val,"/",sep=""))}
if(!dir.exists(fig_dir)){dir.create(fig_dir)}
#-----------------------------------------------------
# trim outliers for dv, otherwise weird-looking qqplot
dv2<-paste(dv,trim_val,sep="_")
d0[,dv2]<-censor(d0[,dv],fraction=trim_val)
w2<-which(is.na(d0[,dv2])) # trim outliers
d0$trimmed<-"no"
d0$trimmed[w2]<-dv2
table(d0$trimmed)
rm(w2,dv2)
  
# save lists of IDs per subset
write.table(d0[,c("subjectid","trimmed")],file=paste(base_dir,"/",v,"_baseline_sMRI_baseline_subjectIDs_",s ,"_trim",trim_val,".txt",sep=""))

# for all samples
summary_all<-summary_vars(d=d0,n_vars=n_vars,c_vars=c_vars)
summary_all$subset<-s
  
# complete cases after trimming dependent variable
summary_trimmed<-summary_vars(d=d0 %>% filter(trimmed=="no"),n_vars=n_vars,c_vars=c_vars)
summary_trimmed$subset<-paste(s, "-", "trimmed",sep="")
  
summary_comb<-rbind(summary_all,summary_trimmed)

# save summary statistics
write.csv(summary_comb,file=paste(base_dir,"/summary_variables_",dv,"_trim",trim_val,".csv",sep=""),row.names = FALSE, quote = TRUE)

# visualize variables and save plots
for (c in vars){
  
  if(class(d0[,c])=="numeric"){
    w2<-which(d0$trimmed=="no")
    p<-ggplot(data=d0) + geom_histogram(aes_string(x=c,fill="trimmed"),bins=50) + 
      geom_vline(aes(xintercept=range(d0[w2,c])[1]),linetype="dotted") +
      geom_vline(aes(xintercept=range(d0[w2,c])[2]),linetype="dotted") +
      scale_fill_brewer(palette="Dark2") +
      theme(legend.position="bottom") +
      labs(subtitle = paste(v %>% gsub("_"," ",.) %>% gsub("ABCD","ABCD ",.),
                            "\n",
                            s," (N = ",NROW(d0),", of which trimmed =",sum(d0$trimmed!="no"),")",sep="")) + 
      NULL
    ggsave(p,file=paste(fig_dir,"/","distribution_",c,".pdf",sep=""))
    rm(p)
    rm(w2)
  }
  
  if(class(d0[,c])=="factor"){
    p<-ggplot(data=d0) + geom_bar(aes_string(x=c,fill="trimmed")) + 
      labs(subtitle = paste( v %>% gsub("_"," ",.) %>% gsub("ABCD","ABCD ",.),
                            "\n",
                            s," (N = ",NROW(d0),", of which trimmed =",sum(d0$trimmed!="no"),")",
                            sep="")) + 
      scale_fill_brewer(palette="Dark2") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme(legend.position="bottom") +
      NULL
    ggsave(p,file=paste(fig_dir,"/","barplot_",c,".pdf",sep=""))
    rm(p)
  }
  
  
}
rm(c)
  
# clean
rm(summary_comb)
rm(summary_trimmed,summary_all)

