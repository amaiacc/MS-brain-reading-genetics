#---------------------------------------------------
# LMM with PRS - r2 per run
#---------------------------------------------------
rm(list=ls())
args <- commandArgs(TRUE)
config_file<-args[1]
# config_file="F:/projects/resources/datasets/ABCD/scripts/genetics/PRS/base_cognitive.config"
# config_file="/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/scripts/genetics/PRS/base_UKB_BIG40.config"
# config_file="F:/projects/resources/datasets/ABCD/scripts/genetics/PRS/base_readingGenLang.config"

# ---------------------------------------------------
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
geno_dir=paste(dir,project,"data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
pheno_dir=paste(dir,project,"/data/working_data/phenotypes/",v,"/",sep="")
input_dir=paste(pheno_dir,s,"/",pheno,"/","baseline","_",atlas,"/","trim",trim_val,"/",sep="") %>% gsub(" ","",.)
prs_dir=paste(dir,project,"/data/working_data/genotyping/prsice/",base_project,"/",sep="")
ukb_dir=paste(dir,"/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")

# define paths for output
out_dir=paste(prs_dir,"/output/",sep="")

# create dirs and set working dir
if(!dir.exists(out_dir)){dir.create(out_dir)}
setwd(prs_dir)

#---------------------------------------------------
cat("---------------------------------\n")
cat("Running subset: ",s,"\n")
cat("PC file: ",pc_file,"\n")
cat("---------------------------------\n")


#---------------------------------------------------
# define ROIs from the brain-behaviour analysis, based on the main analysis, i.e. using subset="all"
rois2follow<-read.table(paste(pheno_dir,"/all/",pheno,"/",atlas,"/summary/tables/","ROImeasures2follow.txt",sep=""),header=FALSE) %>% pull(V1) %>% as.vector()
rois<-rois2follow[grep("total|temp.sup",rois2follow)] # only measures that show some rg with reading
rois2match<-rois %>%  gsub("smri_|area_|thick_|cort.destrieux_","",.) %>% gsub("\\.lh|\\.rh","",.) %>%
  gsub("\\."," ",.) %>%
  unique()
# read base pheno idp info
idps<-read.csv(paste(ukb_dir,"IDPs_summary.csv",sep="")) %>% select(-Units,Type,Cat.,Category_name,downloaded,N.discovery..1,N.replication..1)
idps<-idps %>% filter(Global==1 | ABCD_BrainBehaviour==1)
idps$name2match<-gsub("-"," ",idps$region) %>% gsub("\\+"," and ",.) %>% tolower()

global_idps<-idps %>% filter(Global==1) %>% filter(measure=="AREA") %>% pull(Pheno)

#---------------------------------------------------
# read data
#---------------------------------------------------
# pgs + genetic PCs (for all and EUR only)
pgs<- read.csv(paste(prs_dir,"/QC/","PGS_PCs_",base_project,".csv",sep=""))
#
pgs_cols<-colnames(pgs)[grep("_",colnames(pgs))] %>% unique()
pgs_thresholds<-unique(sapply(strsplit(pgs_cols,"_"),"[[",2))
pgs_thresholds2<-gsub("e.","e-",gsub("^X","",pgs_thresholds))

#---------------------------------------------------
# baseline info/phenotypes
f=paste(input_dir,v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.csv",sep="")
d<-read.csv(f)
d$rel_family_id<-as.factor(d$rel_family_id)

#---------------------------------------------------
# define dependent variables to consider
## brain measures
smri_tot<-colnames(d)[grep("smri*.*destrieux*.*total.lh|smri*.*destrieux*.*mean.lh",colnames(d))]
smri_tot<-smri_tot[grep("sulc|vol|thick",smri_tot,invert=TRUE)]
smri_cols<-colnames(d)[grep("smri*.*destrieux",colnames(d))]

vars2test_b<-smri_cols[grep(paste(rois,collapse="|"),smri_cols)] %>% unique()

## cognitive measures
nihtbx_cols<-colnames(d)[grep("nihtbx*.*uncorrected",colnames(d))]
vars2test_c<-nihtbx_cols[grep("reading",nihtbx_cols)]
#---------------------------------------------------
# define covariates
covs_order<-c("sex","age","high.educ.bl","household.income.bl"
              # "race.4level"
              # ancestry_vars,"hisp","married.bl"
) #
covs_int=c()
# covs_int<-"I(age^2) + sex * age + sex * I(age^2)"

#---------------------------------------------------
# combine
pgs_d<-merge(pgs,d,by.x=c("id"),by.y="subjectid") %>% rename(subjectid=id)
pgs_d$sex<-factor(pgs_d$sex)
pgs_d$rel_family_id<-factor(pgs_d$rel_family_id)

# select columns
d<- pgs_d %>% select("subjectid",
                     "rel_family_id","abcd_site", "mri_info_device.serial.number",
                     contains("PC"),
                     all_of(covs_order),
                     all_of(vars2test_b), all_of(smri_tot),
                     all_of(vars2test_c),
                     all_of(pgs_cols)
                     )

#' #------------------------------------
#' # define number of independent observations
#' #' To correct for multiple testing error, the effective number of independent observations 
#' #' (calculated from a correlation matrix of 8 PGS thresholds x 3 ROIs) was estimated using Matrix 
#' #' Spectral Decomposition (MatSpD) (Nyholt, 2004) before undergoing Bonferroni correction.
#' 
#' # run matSpD to define numbe of indep. obs
#' if(!file.exists("matSpD.out")){
#'   r<-cor(d_pgs[,c(pgs_cols,vars2test)],use="pairwise.complete.obs")
#'   colnames(r)<-rownames(r)<-gsub("\\.","_",colnames(r)) %>% gsub("_","",.)
#'   write.table(r,"correlation.matrix",quote=FALSE,sep="\t")
#'   rm(r) 
#'   source(matSpD_file,local=TRUE, echo=TRUE)
#'   n2correct<-MeffLi
#'   # clean intermediate objects
#'   rm(list=ls(pattern="evals|Old|New|new|old|Message|temp|Warning|Separator|Original|VAR"))
#'   rm(list=ls(pattern="prcoeff|prmax|prindep|row|Meff|vr|ur|svdout|svdsdev|loadings|corr.matrix|coefficients"))
#'   rm(M,L,i,factors,MeffLi)
#' } else {
#'   n2correct<-22 # should read from matSpD.out
#' }
#' 
#' #------------------------------------


# if sumstats exist already, load them, else compute
if (file.exists(paste0(out_dir,"PGS_",base_project,"_",s,"_summary_stats.csv"))){
  pgs_summary<-read.csv(paste0(out_dir,"PGS_",base_project,"_",s,"_summary_stats.csv"))
  pgs_summary$threshold<- pgs_summary$threshold %>% as.character %>% as.numeric %>% as.factor()
  
} else {
  
for(p in c(vars2test_c,vars2test_b,smri_tot) %>% unique()){ # run for each dv
  
  pheno2pcs_cor<-do.call("rbind",lapply(paste("PC",1:10,sep=""),function(pc) {
    z<-cor.test(d[,p],d[,pc])
    zz<-cbind(P=p,PC=pc,estimate=z$estimate,p.value=z$p.value)
    return(zz)
  } )) %>% data.frame()
  
  pcs2keep<-pheno2pcs_cor$PC[which(as.numeric(pheno2pcs_cor$p.value)<0.05)]
  pcs2keep<-paste("PC",1:10,sep="")
  
  # define random depending on dependent variable, brain or cogn
  if (p %in% vars2test_b){
    random1<- "(1|mri_info_device.serial.number)" 
    random2<-"(1|mri_info_device.serial.number) + (1|mri_info_device.serial.number:rel_family_id)"

  } else {
    random1<- "(1|abcd_site)" 
    random2<-"(1|abcd_site) + (1|abcd_site:rel_family_id)"
  }
  
  # define random part
  if(length(grep("unrelated",s))==0){
    random<-random2
  }else{
    ## for unrelated
    random<-random1
  }
  #----------------------
  # define all covariates for the baseline model
  covs2add<-c(covs_order,pcs2keep,covs_int)

  # samples that do not have NA in covariates
  w<-which(apply(d[,covs2add],1,function(x) sum(is.na(x)))==0)

  # z-scores for all the quantiative dependent and independent variables
  d2<-d[w,]  %>% mutate_if(is.numeric,scale)
  d2$rel_family_id<-factor(d2$rel_family_id)
  
  # base formula, including covariates
  # include random effects
  model_formula<-paste(p,"~",paste(c(random,covs2add),collapse=" + "))
  m1<-do.call("lmer",list (as.formula(model_formula),data=d2))

  # define pgs-s to test, which will depend on the specific run...
  if (base_project=="UKB_BIG40" & (p %in% vars2test_b)){
    region<-rois2match[match(p,rois)]
    phenos2run<-idps$Pheno[which(idps$name2match %in% region)]
    pattern2check<-c(phenos2run,global_idps) %>% unique() %>% paste(.,collapse="|",sep="")
    pgs_cols2<-pgs_cols[grep(pattern2check,pgs_cols)]
    rm(region,phenos2run,pattern2check)
  } else {
    pgs_cols2<-pgs_cols
  }
  
  # if p is a brain phenotype, only run the PGS that is best predictor for reading! (not all possibilities)
  if (p %in% vars2test_b){
    if (exists("pgs_summary")){
      tmp<- pgs_summary %>% filter(phenotype %in% vars2test_c)
      maxR2<-max(pgs_summary$R2)
      bestPthr<-tmp %>% filter(R2==maxR2) %>% pull(threshold)
      bestBase<-tmp %>% filter(R2==maxR2) %>% pull(base_pheno)
      bestPGS<-paste(bestBase,bestPthr,sep="_")
      # define PGS to test
      pgs_cols2<-bestPGS
      # clean intermediate
      rm(tmp,maxR2,bestPthr,bestBase,bestPGS)
      
    }
  }
  
  # empty tmp df to be filled with results
  t<-data.frame(phenotype=p,base_model=model_formula,threshold=pgs_cols2,R2=NA,P=NA)
  t[,colnames(summary(m1)$coefficients)]<-NA # # add cols for coefficients

  for (thr in pgs_cols){
      formula_t<-paste(model_formula,"+ ", thr)
      mt<-do.call("lmer",list (as.formula(formula_t),data=d2))
      t$R2[t$threshold==thr]<-(r.squaredGLMM(mt)[1]-r.squaredGLMM(m1)[1])
      t$P[t$threshold==thr]<-anova(m1,mt)$`Pr(>Chisq)`[2]
      t[which(t$threshold==thr),colnames(summary(mt)$coefficients)]<-summary(mt)$coefficients[thr,]
      rm(mt)
    }
  rm(thr,model_formula,m1,pcs2keep,pheno2pcs_cor)
  # if phenotype is brain -> run also correcting for total CSA or mean thickness
  if (p %in% rois){
    # base formula, including covariates
    # include random effects
    if(length(grep("area|vol",p))==1){
      cov_total="smri_area_cort.destrieux_total.lh"
    }
    if(length(grep("thick",p))==1){
      cov_total="smri_thick_cort.destrieux_mean.lh"
    }
    
    model_formula<-paste(p,"~",paste(c(random,covs2add,cov_total),collapse=" + "))
    m2<-do.call("lmer",list (as.formula(model_formula),data=d2))
    
    # empty tmp df to be filled with results
    t2<-data.frame(phenotype=p,base_model=model_formula,threshold=pgs_cols2,R2=NA,P=NA)
    t2[,colnames(summary(m2)$coefficients)]<-NA # # add cols for coefficients
    
    for (thr in pgs_cols2){
      formula_t<-paste(model_formula,"+ ", thr)
      mt<-do.call("lmer",list (as.formula(formula_t),data=d2))
      t2$R2[t$threshold==thr]<-(r.squaredGLMM(mt)[1]-r.squaredGLMM(m2)[1])
      t2$P[t$threshold==thr]<-anova(m2,mt)$`Pr(>Chisq)`[2]
      t2[which(t$threshold==thr),colnames(summary(mt)$coefficients)]<-summary(mt)$coefficients[thr,]
      rm(mt,formula_t)
    }
    rm(thr,model_formula,m2,covs2add,pcs2keep,pheno2pcs_cor)
    t<-rbind(t,t2) 
    rm(t2)
  }

  #
  t$base_pheno<-sapply(strsplit(as.character(t$threshold),"_"),"[[",1) #,sapply(strsplit(as.character(t$threshold),"_"),"[[",2)
  t$threshold<-sapply(strsplit(as.character(t$threshold),"_"),"[[",2)
  t$threshold<-gsub("e.","e-",gsub("^X","",t$threshold))
  t$N<-length(w)
  t$subset<-s
    
  if(exists("pgs_summary")){  pgs_summary<-rbind(pgs_summary,t) } else { pgs_summary<-t }
  rm(t,p,w)
  rm(random)
  
}

rm(d,d2)

pgs_summary$threshold<- pgs_summary$threshold %>% as.character %>% as.numeric %>% as.factor()
pgs_summary$target_pheno<-pgs_summary$phenotype %>%  gsub("nihtbx_|_uncorrected|smri_|_cort.|destrieux","",.)
pgs_summary$adj<-""
pgs_summary$adj[grep("+ smri_vol_cort.destrieux_total",pgs_summary$base_model)]<-"tBV"
pgs_summary$adj[grep("+ smri_area_cort.destrieux_total",pgs_summary$base_model)]<-"tCSA"
pgs_summary$adj[grep("+ smri_thick_cort.destrieux_mean",pgs_summary$base_model)]<-"mCT"


# save output
write.csv(pgs_summary,file=paste0(out_dir,"PGS_",base_project,"_",s,"_summary_stats.csv"),row.names = FALSE)
saveRDS(pgs_d,file=paste0(out_dir,"PGS_",base_project,"_",s,"_data.rds")) 
}
