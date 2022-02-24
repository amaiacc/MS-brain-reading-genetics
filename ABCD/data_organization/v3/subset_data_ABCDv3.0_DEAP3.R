library(dplyr)
library(tidyr)
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
data_dir=paste(dir,"resources/datasets/ABCD/data/primary_data/ABCD3DEAP/results/",sep="")
out_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",sep="")
setwd(data_dir)
# read list of files to be read from the raw txt files
v="ABCDv3_DEAP"
# read data and combine
dat <- readRDS("nda3.0.Rds")

# check structure of dataset
dim(dat)

# define columns of interest
cols<-colnames(dat)
# age sumstats per event (timepoint)
tapply(dat$age,dat$event_name,summary)

# rename some columns, to match v2 colnames
dat<-dat %>% mutate(EUR=genetic_ancestry_factor_european,
                      AFR=genetic_ancestry_factor_african,
                      EAS=genetic_ancestry_factor_east_asian,
                      AMR=genetic_ancestry_factor_american,
                      paired.subjectid=genetic_paired_subjectid_1,
                      PI_HAT=genetic_pi_hat_1,
                      Zygosity=genetic_zygosity_status_1
)


# subset dataset to only baseline
dat1<-subset(dat,eventname=="baseline_year_1_arm_1") # NROW=11878
cols2keep_genetics<-c("subjectid","sex",
                      "EUR","AFR","EAS","AMR",
                      "paired.subjectid",
                      "PI_HAT","Zygosity",
                      cols[grep("genetic_",cols)],
                      cols[grep("^race|ethnicity",cols)],
                      cols[grep("^rel_",cols)]
                      )
dat1[,cols2keep_genetics] %>% 
  write.table(paste(out_dir,v,"_sex_ancestry.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

# structural imaging measures
cols[grep("smri",cols)] %>% strsplit(.,"_") %>% sapply(.,"[[",2) %>% table() # which types of structural MRI measures have been extracted?
cols[grep("smri_area",cols)] %>% strsplit(.,"_") %>% sapply(.,"[[",3) %>% table() # taken e.g. area, which atlases? desikan and destrieux

# define columns that contain rois, for desikan and destrieux atlas
rois_tmp<-cols[grep("desikan|destrieux",cols)]
rois<-rois_tmp[grep("smri_area|smri_thick|smri_vol|smri_sulc",rois_tmp)]
rm(rois_tmp)

# get total brain volume measures, and compute 2/3
dat1[,paste(rois[grep("vol*.*total$",rois)],"2.3",sep="")]<-NA
dat1[,paste(rois[grep("vol*.*total$",rois)],"2.3",sep="")]<-apply(dat1[,rois[grep("vol*.*total$",rois)]],2,function(x) x^(2/3))


# define columns to keep
## could add new ones if needed
cols2keep<-c( "subjectid","eventname",
              "abcd_site",
              "sex","age",
              "rel_family_id",
              "race.4level","race.6level","hisp",
              "AFR","EAS","EUR","AMR", # Proportion of African/East Asian/European/American ancestry
              "paired.subjectid", "PI_HAT","Zygosity", # Twins: subject ID of paired twin that PI_HAT based on, probability of identity by descend (only for twins I think)
              "high.educ.bl","married.bl","household.income.bl", # SES covariates
              cols[grep("nihtbx",cols)], # NIH toolbox related measures, including reading and picture vocabulary
              cols[grep("pea_",cols)],
              cols[grep("neurocog_pc",cols)], # PCs from Thompson et al. 2019
              cols[grep("mri_info",cols)], # info about MRI scanner
              "fsqc_qc",
              rois,
              paste(rois[grep("vol*.*total$",rois)],"2.3",sep="")
)
# make sure that all selected columns are part of the dataset
cols2keep[! cols2keep %in% colnames(dat1) ]

# create subset to save
dat1_s<-dat1 %>% filter(fsqc_qc==1|fsqc_qc=="accept") %>% select(all_of(cols2keep)) # NROW=11265
# save
if(!dir.exists(out_dir)){dir.create(out_dir)}
saveRDS(dat1_s,paste(out_dir,v,"_baseline_sMRI.Rds",sep=""))
