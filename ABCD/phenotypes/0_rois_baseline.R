# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]

config_file="reading.config" # parameters for this run

#---------------------

#---------------------------------------------------
# define libraries
library(dplyr); library(tidyr)
library(data.table)
library(gamm4) #
library(lme4)
library(moments) # skewness and kurtosis
library(MuMIn)
library(psych)
library(lmerTest)
# library(car) # to get vif
#
library(sjPlot)
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
#---------------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
#---------------------
# source config file to define parameters for this run
source(paste(scripts_dir,config_file,sep=""))
#---------------------
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",sep="")
#
working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",sep="")
base_dir=paste(pheno_dir,"baseline","_",atlas,"/","trim",trim_val,"/",sep="")
if(!dir.exists(base_dir)){dir.create(base_dir)}
if(!dir.exists(paste(base_dir,"/tables/",sep=""))){dir.create(paste(base_dir,"/tables/",sep=""))}
if(!dir.exists(paste(base_dir,"/figures/",sep=""))){dir.create(paste(base_dir,"/figures/",sep=""))}
if(!dir.exists(paste(base_dir,"/figures/","/model_diagnosis/",sep=""))){dir.create(paste(base_dir,"/figures/model_diagnosis/",sep=""))}

setwd(working_dir)
#---------------------------------------------------
# read data
dat <- readRDS(paste(working_dir,v,"_baseline_sMRI.Rds",sep=""))
ids<-read.table(file=paste(base_dir,"/../",v,"_baseline_sMRI_baseline_subjectIDs_",s ,"_trim",trim_val,".txt",sep="")) # %>% rename(subjectid=V1,trimmed=V2)
# get PCs, from file defined in config
pcs<-fread(paste(geno_dir,pc_file,sep=""),header=FALSE) ## all samples
colnames(pcs)<-c("FID","subjectid",paste("PC",1:100,sep=""),"pheno")
pcs<-pcs %>% select(-pheno)
## define dataset to be used
d<-merge(dat,ids,by="subjectid")
d<-merge(d,pcs,by="subjectid")
d<-d %>% filter(trimmed=="no")
d$rel_family_id<-as.factor(d$rel_family_id)
rm(dat,ids,pcs)
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
# defines rois to be used in analisis
mri_info<-colnames(d)[grep("mri_info*.*device",colnames(d))]
total<-colnames(d)[grep("total|mean",colnames(d))] # measures are same for desikan
total<-total[grep(atlas,total)]
total_vol<-total[grep("vol*.*total$",total)]
total_vol23<-paste(total_vol,"2.3",sep="")
lh_area<-total[grep("area*.*lh",total)]
lh_thick<-total[grep("thick*.*.lh",total)]


# define rois
rois<-colnames(d)[grep(atlas,colnames(d))]
rois<-rois[grep("thick|area",rois)]
rois<-rois[grep("total",rois,invert=TRUE)]
rois_lh<-rois[grep("lh",rois)]
#---------------------------------------------------
## covariates to consider (for sensitivity)
covs_add<-c(total_vol,cov11,cov12,cov2)

## all variables to consider
vars<-c(covs_def,dv,covs_add,"rel_family_id","abcd_site",
        mri_info,
        rois,
        total,total_vol23
) #ancestry_vars,

#---------------------------------------------------
# define dataset to be used
## complete cases
d<-d %>% select("subjectid",vars)
d<-d[complete.cases(d),]
d<-d %>% mutate_if(is.numeric,scale)
f=paste(base_dir,v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.csv",sep="")
f2=gsub(".csv",".RDS",f)
# write.csv(d,f,row.names=FALSE)
# saveRDS(d, file = f2) # save as RDS as well, as it saves the ordered factor levels for object
rm(f,f2)

#---------------------------------------------------
# define random part
if(length(grep("unrelated",s))==0){
  random<-"(1|mri_info_device.serial.number) + (1|mri_info_device.serial.number:rel_family_id)"
}else{
  ## for unrelated
  random<-"(1|mri_info_device.serial.number)"
}
#---------------------------------------------------
# ## Baseline models
model_b0<-paste(dv,"~",paste(c(random,covs_def),collapse=" + ",sep=" "),sep=" ")

# run baseline model, and do some model diganostics
print(model_b0)
m0<-do.call("lmer",list (as.formula(model_b0),data=d))
m0.gam<-gamm4(as.formula(paste(dv,"~",paste(covs_def,collapse=" + "))),
              random=as.formula(paste("~",random)), data=d,REML=TRUE)

# test<-data.frame(res.m0=residuals(m0),fitted.m0=fitted(m0),
#                  res.m0.mer=residuals(m0.gam$mer),fitted.m0.mer=fitted(m0.gam$mer),
#                  res.m0.gam=residuals(m0.gam$gam),fitted.m0.gam=fitted(m0.gam$gam)
# )
#---------------------------------------------------
# visualize model
## assumptions
m0_plot_diag<-plot_model(m0,type="diag")
m0_plot_diag[5]<-m0_plot_diag[[2]][1]
m0_plot_diag[6]<-m0_plot_diag[[2]][2]
plot_grid(plotlist=m0_plot_diag[c(1,3:6)],nrow=2) %>% 
  ggsave(.,
         file=paste(base_dir,
                    "/figures/model_diagnosis/",dv,"_m0","_jsplot",".png",sep=""),
         width=20,height=6
  )
## slope of coeffs for each predictor vs dv
m0_plot_diag1<-plot_model(m0,type="slope")
ggsave(m0_plot_diag1,
         file=paste(base_dir,
                    "/figures/model_diagnosis/",dv,"_m0","_slopes_jsplot",".png",sep="")
  )
# estimates
m0_plot_estimates<-plot_model(m0,type="est",show.values=TRUE)
ggsave(m0_plot_estimates,
       file=paste(base_dir,
                  "/figures/",dv,"_m0","_estimates_jsplot",".png",sep="")
)
# clean intermediate
rm(list=ls(pattern="m0_plot"))
gc()
#---------------------------------------------------
# Residual plots for baseline model
cat("Residual plot for model0. If it shows structure, it would indiate that there is some deviation from linear form.\n")

png(paste(base_dir,"/figures/model_diagnosis/",dv,"_m0","_residuals",".png",sep=""))
plot(m0) %>% print()
dev.off()

png(paste(base_dir,"/figures/model_diagnosis/",dv,"_m0","_residuals.gam",".png",sep=""))
plot(m0.gam$gam$fitted,m0.gam$gam$residuals)
dev.off()

#' Normality of residuals
qq<-ggplot(data=data.frame(x=residuals(m0)), aes(sample = x)) +  # Create QQplot with ggplot2 package
  stat_qq() +
  stat_qq_line(col = "red")
qq.gam<-ggplot(data=data.frame(x=residuals(m0.gam$gam)), aes(sample = x)) +  # Create QQplot with ggplot2 package
  stat_qq() +
  stat_qq_line(col = "red")
ggsave(qq,file=paste(base_dir,"/figures/model_diagnosis/",dv,"_m0","_qqplot",".png",sep=""))
ggsave(qq,file=paste(base_dir,"/figures/model_diagnosis/",dv,"_m0.gam","_qqplot",".png",sep=""))
rm(qq,qq.gam)

# Linearity
linearity_covs_plots<-lapply(covs_def,function(c){
  if(class(d[,c])=="factor"){
    ggplot(data.frame(c=d[,c],residuals=residuals(m0,type="pearson")),
           aes(x=c,y=residuals)) +
      # geom_violin() + 
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      xlab(c)
  } else {
    ggplot(data.frame(c=d[,c],residuals=residuals(m0,type="pearson")),
           aes(x=c,y=residuals)) +
      geom_point() + geom_smooth(method=lm) +
      xlab(c)
  }
})

png(paste(base_dir,"/figures/model_diagnosis/",dv,"_m0","_residuals_covs",".png",sep=""),width=1000, height=1000)
cowplot::plot_grid(plotlist=linearity_covs_plots) %>% print()
dev.off()
rm(linearity_covs_plots)


#' Sensitivity to data
levp<-ggplot(data.frame(lev=hatvalues(m0),pearson=residuals(m0,type="pearson")), aes(x=lev,y=pearson)) +
      geom_point() + labs(title ="Leverage") +
      NULL
ggsave(levp,file=paste(base_dir,"/figures/model_diagnosis/",dv,"_m0","_leverage",".png",sep=""))
rm(levp)
# Which observations have most leverage?
levId <- which(hatvalues(m0) >= .05)
# if (length(levId)>0) {
#   d[levId,c(dv,covs_def)]
#   summary(d[,c(dv,covs_def)])
#   # check if removing these datapoints changes the model
#   # if random effects are defined
#   m0Lev<-do.call("lmer",list (as.formula(model_b0),data=d[-c(levId),]))
#   m0LevCD <- data.frame(effect=fixef(m0),
#                         change=(fixef(m0Lev) - fixef(m0)),
#                         se=sqrt(diag(vcov(m0)))
#   )
#   rownames(m0LevCD) <- names(fixef(m0Lev))
#   m0LevCD$multiples <- abs(m0LevCD$change / m0LevCD$se)
#   m0LevCD %>% kable(caption="Change in the value of coefficients after removing datapoints with highest leverage.")
#   rm(levId,m0Lev,m0LevCD)
# }
rm(levId)

# save estimates from baseline model
# summary table, parameter table - to get t-values
sum_table<-summary(m0) %>% coef()  %>% data.frame()
# ANOVA table, for all parameters
anova_table<-anova(m0) %>% data.frame() 
# effect size, r2 for the fixed effects(r2m)
effect_table<-data.frame(model="m0",R2m=r.squaredGLMM(m0)[1], R2c=(r.squaredGLMM(m0)[2]))
# save summary tables
write.csv(anova_table,file=paste(base_dir,"/tables/baseline","_","m0","_","anova_table.csv",sep=""),row.names = TRUE)
write.csv(sum_table,file=paste(base_dir,"/tables/baseline","_","m0","_","summary_t_table.csv",sep=""),row.names = TRUE)
write.csv(effect_table,file=paste(base_dir,"/tables/baseline","_","m0","_","effect_r2m_table.csv",sep=""),row.names = FALSE)
# clean
rm(sum_table,anova_table,effect_table)

#---------------------------------------------------

#---------------------------------------------------
# other models
models_table<-data.frame(model_name=c(
                                # m0
                                "m0",
                                "m0b",
                                "m0b2.3",
                                "m0c",
                                "m0d",
                                # m1
                                "m11",
                                "m11b",
                                "m11b2.3",
                                "m11c",
                                "m11d",
                                
                                "m12",
                                "m12b",
                                "m12b2.3",
                                "m12c",
                                "m12d",
                                # m2
                                "m21",
                                "m21b",
                                "m21b2.3",
                                "m21c",
                                "m21d",
                                
                                "m22",
                                "m22b",
                                "m22b2.3",
                                "m22c",
                                "m22d"
                              ),
                              formula=c(
                                # m0
                                model_b0,
                                paste(model_b0, total_vol,sep=" + "),
                                paste(model_b0, total_vol23,sep=" + "),
                                paste(model_b0, lh_area,sep=" + "),
                                paste(model_b0, lh_thick,sep=" + "),
                                # m1
                                paste(model_b0,cov11,sep=" + "),
                                paste(model_b0,cov11,total_vol,sep=" + "),
                                paste(model_b0,cov11,total_vol23,sep=" + "),
                                paste(model_b0,cov11,lh_area,sep=" + "),
                                paste(model_b0,cov11,lh_thick,sep=" + "),
                                
                                paste(model_b0,cov12,sep=" + "),
                                paste(model_b0,cov12,total_vol,sep=" + "),
                                paste(model_b0,cov12,total_vol23,sep=" + "),
                                paste(model_b0,cov12,lh_area,sep=" + "),
                                paste(model_b0,cov12,lh_thick,sep=" + "),
                                # m2
                                paste(model_b0,cov11,cov2,sep=" + "),
                                paste(model_b0,cov11,cov2,total_vol,sep=" + "),
                                paste(model_b0,cov11,cov2,total_vol23,sep=" + "),
                                paste(model_b0,cov11,cov2,lh_area,sep=" + "),
                                paste(model_b0,cov11,cov2,lh_thick,sep=" + "),
                                
                                paste(model_b0,cov12,cov2,sep=" + "),
                                paste(model_b0,cov12,cov2,total_vol,sep=" + "),
                                paste(model_b0,cov12,cov2,total_vol23,sep=" + "),
                                paste(model_b0,cov12,cov2,lh_area,sep=" + "),
                                paste(model_b0,cov12,cov2,lh_thick,sep=" + ")
                                
                              )
)

# save model tables
for (i in 2:NROW(models_table)){
  m<-models_table$model_name[i]
  
  if (file.exists(paste(base_dir,"/tables/baseline","_",m,"_","summary_t_table.csv",sep=""))==FALSE){
  model<-models_table$formula[i] %>% as.character()
  m_t<-do.call("lmer",list (as.formula(model),data=d))
  # summary table, parameter table - to get t-values
  sum_table<-summary(m_t) %>% coef()  %>% data.frame()
  # ANOVA table, for all parameters
  anova_table<-anova(m_t) %>% data.frame() 
  # effect size, r2 for the fixed effects(r2m)
  effect_table<-data.frame(model=m,R2m=r.squaredGLMM(m_t)[1], R2c=(r.squaredGLMM(m_t)[2]))
  # save summary tables
  write.csv(anova_table,file=paste(base_dir,"/tables/baseline","_",m,"_","anova_table.csv",sep=""),row.names = TRUE)
  write.csv(sum_table,file=paste(base_dir,"/tables/baseline","_",m,"_","summary_t_table.csv",sep=""),row.names = TRUE)
  write.csv(effect_table,file=paste(base_dir,"/tables/baseline","_",m,"_","effect_r2m_table.csv",sep=""),row.names = FALSE)
  # clean
  rm(m_t,model,sum_table,anova_table,effect_table)
  }
  rm(m)
}
rm(i)

write.csv(models_table,file=paste(base_dir,"/tables/baseline_models_table.csv",sep=""),row.names = FALSE)

