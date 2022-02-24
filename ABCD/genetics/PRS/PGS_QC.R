#---------------------------------------------------
# QC - visualize distributions of PRS scores, and organize data
#---------------------------------------------------

args <- commandArgs(TRUE)
config_file<-args[1]
# config_file="F:/projects/resources/datasets/ABCD/scripts/genetics/PRS/base_UKB_BIG40.config"
# config_file="/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/scripts/genetics/PRS/base_UKB_BIG40.config"

# reading
library(data.table)
# stats
library(lme4)
library(lmerTest)
library(MuMIn)
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

# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/resources/datasets/"} else {dir="/export/home/acarrion/acarrion/projects/resources/datasets/"}

# source config file to define parameters for this run
source(config_file)
td=paste(project,"_",batch,sep="")

# general paths across scripts
geno_dir1=paste(dir,project,"/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
pgs_dir=paste(dir,project,"/data/working_data/genotyping/prsice/",base_project,"/",sep="")

# define paths for output
out_dir=paste(pgs_dir,"/QC/",sep="")

# create dirs and set working dir
if(!dir.exists(out_dir)){dir.create(out_dir)}
setwd(pgs_dir)

#---------------------------------------------------
# read data
#---------------------------------------------------

# read genetic PCs to include as covariates
## for each subset
pcs_all<-fread(paste(geno_dir1,td,".pop_strat.pca.evec",sep=""),header=FALSE) ## all samples
colnames(pcs_all)<-c("FID","IID",paste("all","PC",1:100,sep=""),"pheno")
pcs_all$id<-paste(pcs_all$FID,pcs_all$IID,sep="_")

## eur only, after excluding outliers
pcs_eur<-fread(paste(geno_dir1,td,"_EUR_6sd",".pop_strat_outliers.pca.evec",sep=""),header=FALSE) ## all samples
colnames(pcs_eur)<-c("FID","IID",paste("eur","PC",1:100,sep=""),"pheno")
pcs_eur$id<-paste(pcs_eur$FID,pcs_eur$IID,sep="_")
#---------------------------------------------------
# individual PGS 
# read individual PGS for each base dataset
base_ds<-list.files(pgs_dir,pattern="prsice") %>% gsub(".prsice","",.) %>% gsub(paste(td,"_",sep=""),"",.)
pgs0<-data.frame(FID=NA,IID=NA,batch=NA) # empty df to fill in 
for(bd in base_ds){
  pgs_t<-pgs0
  for(td in td){
    tmp<-fread(paste0(pgs_dir,td,"_",bd,".all.score"),header=TRUE)
    tmp$batch<-td
    # merge
    pgs_t<-merge(pgs_t,tmp,all=TRUE)
    # clean
    rm(tmp)
  }
  # rename columns
  tmp_name<-gsub("WR_","WR.",bd) %>% strsplit(.,"_") %>% sapply("[[",1)
  colnames(pgs_t)[4:(NCOL(pgs_t))]<-paste(tmp_name,colnames(pgs_t)[4:(NCOL(pgs_t))],sep="_")
  colnames(pgs_t)<-gsub("-",".",colnames(pgs_t))
  # merge across base datasets
  if(exists("pgs")){
    pgs<-merge(pgs,pgs_t,all=TRUE)
  } else {
    pgs<-pgs_t
  }
  # clean intermediate
  rm(tmp_name,pgs_t)
}
rm(pgs0,bd)

# edit id
pgs<-subset(pgs,!is.na(IID))
# pgs$id<-pgs$IID %>% gsub("`","",.) %>% strsplit(.,"NDAR_")
pgs$id<-pgs$IID %>% gsub("'","",.) %>% strsplit(.,"NDAR") %>% sapply("[[",2) %>% paste("NDAR_",.,sep="") %>% gsub("__","_",.)

# edit colnames
# colnames(pgs)<-gsub("EA3","EA",colnames(pgs))
pgs_cols<-colnames(pgs)[grep("_",colnames(pgs))] %>% unique()
pgs_thresholds<-unique(sapply(strsplit(pgs_cols,"_"),"[[",2))
pgs_thresholds2<-gsub("e.","e-",gsub("^X","",pgs_thresholds))
# rm(pgs_cols)

# combine pgs and genetic PCs
pcs<-merge(pcs_all,pcs_eur,by=c("id","FID","IID"),all=TRUE)
pgs<-merge(pgs,pcs[,c("IID",paste0("all","PC",1:10),paste0("eur","PC",1:10))],by.x="id",by.y="IID")
rm(pcs)

# normalize PGSes
pgs_z<- pgs %>% mutate_if(is.numeric,scale)


#---------------------------------------------------
# save combined PGS result
write.csv(pgs,paste(out_dir,"PGS_PCs_",base_project,".csv",sep=""),row.names=FALSE)
write.csv(pgs_z,paste(out_dir,"zPGS_zPCs_",base_project,".csv",sep=""),row.names=FALSE)

#---------------------------------------------------
# check distribution of polygenic scores
pgs_long<- pgs %>% gather (threshold,PGS,contains("_")) %>%
  mutate(base_pheno=sapply(strsplit(threshold,"_"),"[[",1),
         threshold=gsub(unique(base_pheno) %>% paste(.,"_",sep="",collapse="|"),"",threshold) %>%
           gsub("^X","",.) %>% gsub("e.","e-",.) %>% as.numeric %>% as.factor()
  )
  
pgsz_long<- pgs_z %>% gather (threshold,PGS,contains("_")) %>%
  mutate(base_pheno=sapply(strsplit(threshold,"_"),"[[",1),
         threshold=gsub(unique(base_pheno) %>% paste(.,"_",sep="",collapse="|"),"",threshold) %>%
           gsub("^X","",.) %>% gsub("e.","e-",.) %>% as.numeric %>% as.factor()
  )

pgs_dist <- ggplot(pgs_long, aes(y=threshold,x=PGS,fill=threshold) ) + geom_density_ridges2(scale=4) +
  facet_grid(.~base_pheno,scales="free") +
  scale_fill_brewer(palette = 2) +
  NULL

pgsz_dist <- ggplot(pgsz_long, aes(y=threshold,x=PGS,fill=threshold) ) + geom_density_ridges2(scale=4) +
  facet_grid(.~base_pheno,scales="free") +
  scale_fill_brewer(palette = 2) +
  labs(x="z PGS") +
  NULL

pgs_legend<-get_legend(pgs_dist)

combined_dist<-plot_grid(pgs_dist + theme(legend.position="none"),
                         pgsz_dist + theme(legend.position="none"),
                         pgs_legend,
                         nrow=1,
                         rel_widths = c(1,1,0.2))

# ggsave(pgs_dist,file=paste(out_dir,"PGS_",base_project,"_distributions.png",sep=""),width=7,height=5)
# ggsave(pgsz_dist,file=paste(out_dir,"zPGS_",base_project,"_distributions.png",sep=""),width=7,height=5)
# ggsave(combined_dist,file=paste(out_dir,"PGSzPGS_",base_project,"_distributions.png",sep=""),width=12,height=5)

#---------------------------------------------------
# relationship of PGS to PCs
# pgs_longer <- pgs_long %>% gather (PC,PCv,contains("PC"))
# ggplot(subset(pgs_longer,PC=="allPC1"|PC=="allPC2"), aes(y=PCv,x=PGS,color=threshold) ) +
#   geom_point() +
#   facet_wrap(base_pheno~PC,scales="free") +
#   scale_color_brewer(palette = "Dark2") +
#   NULL
# rm(pgs_longer)
#---------------------------------------------------
