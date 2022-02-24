# clean workspace
rm(list=ls())

#---------------------
args <- commandArgs(TRUE)
config_file<-args[1]
i=args[2]
#
config_file="reading.config" # parameters for this run
i=1 # args[2]
#---------------------------------------------------
# define libraries
library(dplyr); library(tidyr)
library(data.table)
library(psych)
#
library(scales)
library(ggplot2)
# if (!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("caijun/ggcorrplot2")
# library(ggcorrplot2)
library(ggcorrplot)
library(corrplot)
library(cowplot)

#
corrtheme <- theme(
  axis.title.x = element_text(angle = 0, color = 'grey20',face="bold", size=12),
  axis.title.y = element_text(angle = 90, color = 'grey20',face="bold", size=12),
  legend.title=element_text(face="bold", size=12),
  plot.title = element_text(angle = 0, face="bold", size=14),
  axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 10),
  axis.text.y = element_text(angle=0, hjust = 1,vjust=1,size = 10)
) 
#---------------------------------------------------
# define directories
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"} #"/bcbl/home/home_a-f/acarrion/acarrion/projects/"
scripts_dir=paste(dir,"/resources/datasets/ABCD/scripts/phenotypes/",sep="")
source(paste0(scripts_dir,"general_functions.R"))
#---------------------
# source config file to define parameters for this run
source(paste(scripts_dir,config_file,sep=""))
# source(config_file)

#---------------------
geno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/preImputation_QC/intermediate_files/sampleQC/PopStr/",sep="")
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/",v,"/",s,"/",pheno,"/",sep="")
input_dir=paste(pheno_dir,"baseline","_",atlas,"/",sep="") %>% gsub(" ","",.)
working_dir=paste(pheno_dir,atlas,sep="")
if(!dir.exists(working_dir)){dir.create(working_dir)}
setwd(working_dir)
#---------------------------------------------------

#---------------------------------------------------
# read data
f=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis.csv",sep="")
d<-read.csv(f)
d$rel_family_id <- as.factor(d$rel_family_id)
colnames(d) <- colnames(d) %>% gsub("smri_|cort.destrieux_","",.) %>% gsub("\\.|_","_",.)
# residuals
f2=paste(input_dir,"/trim",trim_val,"/",v,"_baseline_sMRI","_",pheno,"_","trim",trim_val,"_","scaled4analysis_ROIs2follow_residuals.csv",sep="")
d_res<-read.csv(f2)
colnames(d_res) <- colnames(d_res) %>% gsub("smri_|cort.destrieux_","",.) %>% gsub("\\.|_","_",.)

rm(f,f2)
#---------------------------------------------------
# define rois
rois2follow <- read.table(paste(pheno_dir,atlas,"summary","tables","ROImeasures2follow.txt",sep="/")) %>%
  mutate(V1 = V1%>% gsub("smri_|cort.destrieux_","",.) %>% gsub("\\.|_","_",.) )

rois_rh<-gsub("_lh","_rh",rois2follow$V1) # only the ones that show assoc. with LH counterpart
rois_lh<- rois2follow$V1 %>% as.character()
rois<-c(rois_lh,rois_rh)
# covs
covs<-colnames(d)[grep("total.|mean.",colnames(d))]
covs<-covs[grep("area|thick",covs)]
#---------------------------------------------------

#---------------------------------------------------
# correlations across rois2follow + visualize
#---------------------------------------------------
dv2="area_total_lh"
dv1="area_g_parietal_sup_lh"
ggplot(data=d,aes_string(x=dv1,y=dv2,color="high_educ_bl")) + 
  # stat_bin_hex(bins = 100) +
  # geom_point(alpha=0.5) +
  stat_density2d() +
  geom_smooth(method="lm", formula=y~x) +
  # facet_grid(high_educ_bl~sex) +
  NULL



#---------------------------------------------------
# RAW measures
#---------------------------------------------------
d1<-d[,rois]
colnames(d1) <- colnames(d1) 
corROIS<-corr.test(d1,minlength =100,adjust="bonferroni")

## to long format
corROIs_long<-do.call("cbind",
  lapply(c("r","se","p"), function(m){
    x <- corROIS[m] %>% as.data.frame()
    x2 <- x %>%
      mutate(measure1=rownames(x)) %>%
      pivot_longer(cols=colnames(x),names_to="measure2",values_to=m)
    return(x2)
    } )
  ) %>% subset(select=which(!duplicated(names(.)))) %>% 
  mutate(
    measure2=gsub("^r.","",measure2),
    n=corROIS$n, adj=corROIS$adjust,type="raw measures") %>%
  # remove rows of phenotypes r with itself
  filter(measure1!=measure2) %>%
  distinct()

tmp1<- corROIs_long %>% rename(m1=measure1, m2=measure2)
tmp2<- corROIs_long %>% rename(m1=measure2, m2=measure1)
tmp<-rbind(tmp1,tmp2)
tmp[which(duplicated(tmp)),] %>% View()
## plot
c<-corROIS$r
p<-corROIS$p

# # left hemi: upper triangle, 
c2<-c[rois_lh,rois_lh]
p2<-p[rois_lh,rois_lh]
# right hemi: lower triangle
c2[lower.tri(c2)]<-c[rois_rh,rois_rh] [c[rois_rh,rois_rh] %>% lower.tri()]
p2[lower.tri(p2)]<-p[rois_rh,rois_rh] [p[rois_rh,rois_rh] %>% lower.tri()]

# diagonal: left and right
diag(c2)<-c[rois_lh,rois_rh]%>% diag()
diag(p2)<-p[rois_lh,rois_rh] %>% diag()

# change names
names<-colnames(c2) %>% gsub(".lh|.rh","",.) %>% gsub("\\.|_"," ",.) %>%
  sapply(simpleCap) %>%
  gsub("Area","CSA",. ) %>% gsub("Thick","CT",.) #
colnames(c2)<-rownames(c2)<-names
colnames(p2)<-rownames(p2)<-names
rm(names)

# ROIs_corplot <- corrplot(corr=c2, p.mat=p2,
#               method="square",
#               order="hclust",
#               addrect=2,
#               # style
#               insig="blank",
#               addCoef.col=TRUE,
#               # addCoefasPercent=TRUE,
#               # type="upper",
#               diag=TRUE,
#               tl.col="black",
#               tl.srt=40,
#               tl.cex=0.7, # control size of the row and colnames
#               number.cex=0.8,
#               number.digits=3,
#               mar=c(0,0,1,0),
#               title="Unadjusted measures"
#               )

corrplot_unadj <- ggcorrplot::ggcorrplot(c2,
           p.mat=p2,
           method="square",
           insig="blank",
           hc.order=TRUE,
           lab=TRUE,
           type="full",
           show.diag=TRUE,
           digits=2
           # colors = c("#6D9EC1", "white", "#E46726"),
           ) + 
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"),
                       midpoint = 0,
                       limits=c(-1,1)) +
  corrtheme +
  labs(x = 'LEFT', y = 'RIGHT',
       fill="Corr. Coef.",
       title = "Correlations",
       subtitle="Raw measures\n") + 
  NULL

# subset matrices, and rename
d2<-corROIS$ci
d2$phenos<-rownames(d2)
d2$n<-corROIS$n
d2 <- d2 %>% mutate(
  region1 = strsplit(phenos,"-",.) %>% sapply("[[",1),
  region2 = strsplit(phenos,"-",.) %>% sapply("[[",2)
) %>% select(-phenos) %>% select(region1,region2,n,p,r,upper,lower)
rownames(d2)<-1:NROW(d2)

# clean intermediate
rm(c2,p2)
rm(corROIS,d1)

#---------------------------------------------------
# Partial correlations: adjusting for total CSA or mean CT
#---------------------------------------------------

# compute partial correlations, per hemisphere and measure type
all_pcorsROIr<-lapply(covs, function(cov){
  measure<-strsplit(cov,"_") %>% sapply(.,"[[",1)
  hemi<-strsplit(cov,"mean.|total.") %>% sapply(.,"[[",2)
  rois2<-rois[grep(paste(measure,hemi,sep="*.*"),rois)]
  rois2<-rois2[grep("total|mean",rois2,invert=TRUE)]
  # compute
  r<-partial.r(d[,c(rois,covs)],rois2,cov)
  corrmatp<-corr.p(r,n=NROW(d),minlength =100,adjust="bonferroni")
  l<-list(r=r,p=corrmatp$p,corrtable=corrmatp$ci)
  return(l)
})

names(all_pcorsROIr)<-covs

# create table
d3<-lapply(1:length(all_pcorsROIr),function(i) all_pcorsROIr[[i]]$corrtable) %>% do.call("rbind",.)
d3$phenos<-rownames(d3)
d3$n <- NROW(d)
d3 <- d3 %>% mutate(
  region1 = strsplit(phenos,"-",.) %>% sapply("[[",1),
  region2 = strsplit(phenos,"-",.) %>% sapply("[[",2)
) %>% select(-phenos) %>% select(region1,region2,n,p,r,upper,lower)
rownames(d3)<-1:NROW(d3)


## build matrix for partial correlations...
# create matrix and fill in with partial correlations
c3<-c[rois_lh,rois_lh] # get matrix structure from total matrix
c3[1:NROW(c3),1:NCOL(c3)]<-0 # blank all

p3<-c[rois_lh,rois_lh] # get matrix structure from total matrix
p3[1:NROW(c3),1:NCOL(c3)]<-1 # blank all


# fill in from the partial correlations
## left hemisphere -> upper diagonal
n<-colnames(all_pcorsROIr$thick_mean_lh$r)
c3[n,n][upper.tri(c3[n,n])]<- all_pcorsROIr$thick_mean_lh$r[upper.tri(all_pcorsROIr$thick_mean_lh$r)]
p3[n,n][upper.tri(p3[n,n])]<- all_pcorsROIr$thick_mean_lh$p[upper.tri(all_pcorsROIr$thick_mean_lh$p)]

rm(n)
n<-colnames(all_pcorsROIr$area_total_lh$r)
c3[n,n][upper.tri(c3[n,n])]<- all_pcorsROIr$area_total_lh$r[upper.tri(all_pcorsROIr$area_total_lh$r)]
p3[n,n][upper.tri(p3[n,n])]<- all_pcorsROIr$area_total_lh$p[upper.tri(all_pcorsROIr$area_total_lh$p)]

rm(n)
## right hemisphere --> lower diagonal
n<-colnames(all_pcorsROIr$thick_mean_rh$r) %>% gsub("rh","lh",.)
c3[n,n][lower.tri(c3[n,n])]<- all_pcorsROIr$thick_mean_rh$r[lower.tri(all_pcorsROIr$thick_mean_rh$r)]
p3[n,n][lower.tri(p3[n,n])]<- all_pcorsROIr$thick_mean_rh$p[lower.tri(all_pcorsROIr$thick_mean_rh$p)]

rm(n)
n<-colnames(all_pcorsROIr$area_total_rh$r) %>% gsub("rh","lh",.)
c3[n,n][lower.tri(c3[n,n])]<- all_pcorsROIr$area_total_rh$r[lower.tri(all_pcorsROIr$area_total_rh$r)]
p3[n,n][lower.tri(p3[n,n])]<- all_pcorsROIr$area_total_rh$p[lower.tri(all_pcorsROIr$area_total_rh$p)]

rm(n)
##
# rename
# change names
names<-colnames(c3) %>% gsub(".lh|.rh","",.) %>% gsub("\\.|_"," ",.) %>%
  sapply(simpleCap) %>%
  gsub("Area","CSA",. ) %>% gsub("Thick","CT",.) #
colnames(c3)<-rownames(c3)<-names
colnames(p3)<-rownames(p3)<-names
rm(names)

# order to match plot of unadjusted...
names3<-levels(corrplot_unadj$data$Var1) %>% as.character()
c3<-c3[names3,names3]
p3<-p3[names3,names3]
rm(names3)


corrplot_adj <- ggcorrplot::ggcorrplot(c3,
                                         p.mat=p3,
                                         method="square",
                                         insig="blank",
                                         hc.order=FALSE,
                                         lab=TRUE,
                                         type="full",
                                         show.diag=TRUE,
                                         digits=2
                                         # colors = c("#6D9EC1", "white", "#E46726"),
) + 
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"),
                       midpoint = 0,
                       limits=c(-1,1)) +
  corrtheme +
  labs(x = 'LEFT', y = 'RIGHT',
       fill="Corr. Coef.",
       title = "Partial correlations",
       subtitle=c("CSA adjustement: total CSA\nCT adjustement: mean CT")) + 
  NULL

# clean intermediate
rm(c3,p3)
rm(all_pcorsROIr)
rm(c,p)


#---------------------------------------------------
# Residuals
#---------------------------------------------------
d1<-d_res[,paste0("residuals_",rois)]
colnames(d1) <- colnames(d1) %>% gsub("residuals_","",.)
corROIS<-corr.test(d1,minlength =100,adjust="bonferroni")
#
c<-corROIS$r
p<-corROIS$p

# # left hemi: upper triangle, 
c4<-c[rois_lh,rois_lh]
p4<-p[rois_lh,rois_lh]
# right hemi: lower triangle
c4[lower.tri(c4)]<-c[rois_rh,rois_rh] [c[rois_rh,rois_rh] %>% lower.tri()]
p4[lower.tri(p4)]<-p[rois_rh,rois_rh] [p[rois_rh,rois_rh] %>% lower.tri()]

# diagonal: left and right
diag(c4)<-c[rois_lh,rois_rh]%>% diag()
diag(p4)<-p[rois_lh,rois_rh] %>% diag()

# change names
names<-colnames(c4) %>% gsub(".lh|.rh","",.) %>% gsub("\\.|_"," ",.) %>%
  sapply(simpleCap) %>%
  gsub("Area","CSA",. ) %>% gsub("Thick","CT",.) #
colnames(c4)<-rownames(c4)<-names
colnames(p4)<-rownames(p4)<-names
rm(names)

order<-levels(corrplot_unadj$data$Var1) %>% as.character()

c4<-c4[order,order]
p4<-p4[order,order]

corrplot_res <- ggcorrplot::ggcorrplot(c4,
                                         p.mat=p4,
                                         method="square",
                                         insig="blank",
                                         hc.order=FALSE,
                                         lab=TRUE,
                                         type="full",
                                         show.diag=TRUE,
                                         digits=2
                                         # colors = c("#6D9EC1", "white", "#E46726"),
) + 
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"),
                       midpoint = 0,
                       limits=c(-1,1)) +
  corrtheme +
  labs(x = 'LEFT', y = 'RIGHT',
       fill="Corr. Coef.",
       title = "Correlations residuals",
       subtitle="Covariates\n") + 
  NULL

# subset matrices, and rename
d4<-corROIS$ci
d4$phenos<-rownames(d4)
d4$n<-corROIS$n
d4 <- d4 %>% mutate(
  region1 = strsplit(phenos,"-",.) %>% sapply("[[",1),
  region2 = strsplit(phenos,"-",.) %>% sapply("[[",2)
) %>% select(-phenos) %>% select(region1,region2,n,p,r,upper,lower)
rownames(d4)<-1:NROW(d4)

# clean intermediate
rm(c4,p4)
rm(corROIS,d1)
rm(c,p)

#---------------------------------------------------
# Residuals adjusted for global measures
#---------------------------------------------------
d1<-d_res[,paste0("residuals_globalAdj_",rois)]
colnames(d1) <- colnames(d1) %>% gsub("residuals_globalAdj_","",.)
corROIS<-corr.test(d1,minlength =100,adjust="bonferroni")
#
c<-corROIS$r
p<-corROIS$p

# # left hemi: upper triangle, 
c5<-c[rois_lh,rois_lh]
p5<-p[rois_lh,rois_lh]
# right hemi: lower triangle
c5[lower.tri(c5)]<-c[rois_rh,rois_rh] [c[rois_rh,rois_rh] %>% lower.tri()]
p5[lower.tri(p5)]<-p[rois_rh,rois_rh] [p[rois_rh,rois_rh] %>% lower.tri()]

# diagonal: left and right
diag(c5)<-c[rois_lh,rois_rh]%>% diag()
diag(p5)<-p[rois_lh,rois_rh] %>% diag()

# change names
names<-colnames(c5) %>% gsub(".lh|.rh","",.) %>% gsub("\\.|_"," ",.) %>%
  sapply(simpleCap) %>%
  gsub("Area","CSA",. ) %>% gsub("Thick","CT",.) #
colnames(c5)<-rownames(c5)<-names
colnames(p5)<-rownames(p5)<-names
rm(names)

order<-levels(corrplot_unadj$data$Var1) %>% as.character()

c5<-c5[order,order]
p5<-p5[order,order]

corrplot_resAdjGlobal <- ggcorrplot::ggcorrplot(c5,
                                                p.mat=p5,
                                                method="square",
                                                insig="blank",
                                                hc.order=FALSE,
                                                lab=TRUE,
                                                type="full",
                                                show.diag=TRUE,
                                                digits=2
                                                # colors = c("#6D9EC1", "white", "#E46726"),
) + 
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"),
                       midpoint = 0,
                       limits=c(-1,1)) +
  corrtheme +
  labs(x = 'LEFT', y = 'RIGHT',
       fill="Corr. Coef.",
       title = "Correlations residuals",
       subtitle="CSA covariate: total CSA\nCT covariate: mean CT\n") + 
  NULL

# subset matrices, and rename
d5<-corROIS$ci
d5$phenos<-rownames(d5)
d5$n<-corROIS$n
d5 <- d5 %>% mutate(
  region1 = strsplit(phenos,"-",.) %>% sapply("[[",1),
  region2 = strsplit(phenos,"-",.) %>% sapply("[[",2)
) %>% select(-phenos) %>% select(region1,region2,n,p,r,upper,lower)
rownames(d5)<-1:NROW(d5)

# clean intermediate
rm(c5,p5)
rm(corROIS,d1)
rm(c,p)

#---------------------------------------------------

out_dir=paste0(working_dir,"/summary/figures/")

plot_grid(corrplot_unadj + theme(legend.position = "none"),
          corrplot_adj + theme(legend.position = "none"),
          corrplot_res + theme(legend.position = "none"),
          corrplot_resAdjGlobal + theme(legend.position = "none"),
          align="hv",
          nrow=1) %>%
  ggsave(file=paste0(out_dir,"corrplot_all_comparison.pdf"))


corrplot_unadj %>%
  ggsave(file=paste0(out_dir,"corrplot_all_unadj.pdf"))

corrplot_adj %>%
  ggsave(file=paste0(out_dir,"corrplot_all_adj.pdf"))

corrplot_res %>%
  ggsave(file=paste0(out_dir,"corrplot_all_res.pdf"))

corrplot_resAdjGlobal %>%
  ggsave(file=paste0(out_dir,"corrplot_all_resAdjGlobal.pdf"))

leg<-get_legend(corrplot_res) %>% plot_grid()
plot_grid(
  plot_grid(corrplot_res + theme(legend.position = "none") + labs(title="",subtitle=""),
          corrplot_resAdjGlobal + theme(legend.position = "none")  + labs(title=""),
          align="hv",
          nrow=1,
          labels=c("A","B")),
  NULL,
  leg,
  rel_widths = c(1,0.05,0.1),nrow=1
  ) %>%
  ggsave(file=paste0(out_dir,"corrplot_res_resAdjGlobal.pdf"),width=13,height=6.2)
