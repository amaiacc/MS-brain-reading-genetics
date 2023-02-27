# quick and dirty comparison of heritability runs
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(MetBrewer)
#---------------------------------------------------
options(stringsAsFactors = FALSE)
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
base_dir=paste0(dir,"resources/datasets/ABCD/data/working_data/genotyping/GCTA/")
setwd(base_dir)
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
#----------------------------------------------------------------------
# read ggseg images and color table, for consistency
pheno_dir=paste(dir,"resources/datasets/ABCD/data/working_data/phenotypes/ABCDv3_DEAP/all/reading/destrieux/summary/",sep="")
load(paste(pheno_dir,"ggseg_brain_ROIs.Rdata",sep=""))
braincolor_codes<- braincolor_codes %>% mutate(region_name=Region  %>% gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% toupper() )
rm(list=ls(pattern="gg_|summary_|title|legend"))
#----------------------------------------------------------------------
# define pattern for this run
root="ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned_rm05_adj"
file=paste0("hsq_h2_summary_",root,".table")
# read output table for each run, and flag them
main<-read.table(paste0("reml_v1/",file)) %>%
  mutate(run="main",subset="GRM",covs="sex,age,PC1:PC10",regression="LMM")
main_lm<-read.table(paste0("reml_v1_lm/",file)) %>%
  mutate(run="main_lm",subset="GRM",covs="sex,age,PC1:PC10",regression="LM")
pheno_subset<-read.table(paste0("reml_v2/",file)) %>%
  mutate(run="pheno_subset",subset="pheno",covs="sex,age,PC1:PC10",regression="LMM")
ses_covs<-read.table(paste0("reml_v3/",file)) %>%
  mutate(run="ses_covs",subset="GRM",covs="sex,age,income,educ,PC1:PC10",regression="LMM")
ses_covs_lm<-read.table(paste0("reml_v3_lm/",file)) %>%
  mutate(run="ses_covs_lm",subset="GRM",covs="sex,age,income,educ,PC1:PC10",regression="LM")
# combine all
h2<-rbind(main,main_lm,pheno_subset,ses_covs,ses_covs_lm)
# rm intermediate
rm(main,main_lm,pheno_subset,ses_covs,ses_covs_lm)
#----------------------------------------------------------------------
# format db to get all columns of interest
colnames(h2)<-c("file","n","h2","se","pval","run","subset","covs","regression")
h2$file<-gsub(":n","",h2$file)
h2$p1<-gsub(".hsq","",h2$file) %>% gsub("_residuals|_covs","",.) 
h2$estimate<-h2$h2
h2$stat<-"h2"
h2$sample<-strsplit(root,"QCed_") %>% sapply("[[",2)
h2$adj<-"-"
h2$adj[grep("adjGlobal",h2$file)]<-"adjGlobal"
# define categories for the phenotypes
h2$category<-NA
h2$category[grep("smri*.*total",h2$p1)]<-"Brain global"
h2$category[grep("smri",h2$p1)[grep("total",h2$p1[grep("smri",h2$p1)],invert=TRUE)]]<-"Brain ROI"
h2$category[grep("smri",h2$p1,invert=TRUE)]<-"Cognitive"
# for brain regions: measure type, hemisphere and region name
h2$measure<-h2$region<-h2$hemisphere<-NA
w<-which(h2$category!="Cognitive")
h2$hemisphere[w][grep("lh|LH",h2$p1[w])]<-"L"
h2$hemisphere[w][grep("rh|RH",h2$p1[w])]<-"R"
h2$region[w]<-gsub(".lh|.rh","",h2$p1[w])
h2$measure[w][grep("area",h2$p1[w])]<-"AREA"
h2$measure[w][grep("thick",h2$p1[w])]<-"THICKNESS"
rm(w)
#
h2 %>% mutate(region_name=region,estimate=h2)
h2$sig<-""
h2$sig[h2$pval<0.05]<-"*"
h2$sig[h2$pval<0.05/NROW(h2)]<-"***"
# define categories for the phenotypes
h2$category<-NA
h2$category[grep("smri*.*total",h2$p1)]<-"Brain global"
h2$category[grep("smri",h2$p1)[grep("total",h2$p1[grep("smri",h2$p1)],invert=TRUE)]]<-"Brain ROI"
h2$category[grep("smri",h2$p1,invert=TRUE)]<-"Cognitive"

# clean phenotype names
h2$pheno1<-h2$p1 %>% gsub("nihtbx_|pea_|_tss|_uncorrected|smri_|_cort.|destrieux|_adjGlobal","",.)

# clean region name
h2$region<-gsub(".rh|.lh","",h2$pheno1)
h2$region_name <- h2$region %>% gsub("smri_thick_cort.|smri_area_cort.|destrieux_","",.) %>% 
  gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% gsub("AREA_|area_|THICK_|thick_","",.) %>% toupper()
h2$measure[grep("area",h2$pheno1)]<-"AREA"
h2$measure[grep("thick",h2$pheno1)]<-"THICKNESS"

h2<-merge(h2,braincolor_codes%>% select(-hemisphere),by=c("region_name","measure"),all.x=TRUE) 
h2<- h2 %>% mutate(run_name=paste0(
  "subset = ",subset,"; regression = ",regression)
  )
h2 <- h2 %>% mutate(name=if_else(str_detect(p1,"WR|reading"),"Word Reading",
                                 if_else(str_detect(p1,"IQ|fluid"),"Fluid",
                                         if_else(str_detect(p1,"wiscv"),"Matrix reasoning",
                                                 if_else(str_detect(p1,"CP"),"Cognitive Performance",
                                                         if_else(str_detect(p1,"cryst"),"Crystalized",
                                                                 if_else(str_detect(p1,"picvocab"),"Vocabulary",
                                                                         if_else(str_detect(p1,"Dyslexia"),"Dyslexia",
                                                                                 if_else(str_detect(p1,"EA"),"Educational Attainement","-"
                                                                                 ))))))))) %>%
  mutate(name2=if_else(str_detect(p1,"WR|reading|Dyslexia"),"Reading",
                       if_else(str_detect(p1,"IQ|CP|cryst|fluid|wisc"),"Intelligence",
                               if_else(str_detect(p1,"picvocab"),"Vocabulary",
                                       if_else(str_detect(p1,"EA"),"Education","-"
                                       ))))) 
# filter phenotypes that are not of interest
h2 <- h2 %>% filter(region_name!="CRYST"&region_name!="MEAN"&
                      ((hemisphere=="L")|category=="Cognitive"))
#----------------------------------------------------------------------
## BRAIN ##
h2_brain<-subset(h2,category!="Cognitive") %>%
  filter(region_name %in% braincolor_codes$region_name) %>% 
  mutate(region_name=factor(region_name,levels=unique(braincolor_codes$region_name))) %>% arrange(region_name) %>%
  mutate(region_name2=factor(region_name %>% as.character() %>% tolower() %>% sapply(.,simpleCap),
                             levels=unique(braincolor_codes$region_name)%>% tolower() %>% sapply(.,simpleCap) )
  )

## COGNITIVE ##
h2_cogn <- subset(h2,category=="Cognitive") %>% distinct() %>%
  filter(name!="-") %>%
  mutate(name2=factor(name2,levels=c("Reading","Vocabulary","Intelligence","Education"))) %>% 
  mutate(name=factor(name,levels=c("Word Reading","Dyslexia","Vocabulary",
                                   "Crystalized",
                                   "Fluid","Matrix reasoning",
                                   "Cognitive Performance",
                                   "Educational Attainement"))) %>% 
  arrange(name2,name)
h2_cogn$name<-factor(h2_cogn$name)


#----------------------------------------------------------------------
# PLOT BRAIN
# compare each run to main analysis
h2_brain_comparison <- ggplot(h2_brain,
         aes(x=region_name2,y=h2,width=0.7,
             # color=region_name,fill=region_name,
             color=covs,fill=covs,
             shape=run_name
             # alpha=sig
             )) + 
  # annotate area and thickness, and add line in between
  annotate("text", x=7-0.25,y=1,label="CSA",size=3) +
  annotate("text", x=6+0.25,y=1,label="CT",size=3) +
  geom_vline(xintercept = seq(1.5,9.5,by=1), color="grey",linetype="solid") +
  geom_vline(xintercept = 7-0.5,color="black",alpha=0.8,linetype="solid",size=1.2) +
  #
  geom_hline(yintercept = 0,color="black",alpha=0.3,linetype="dashed") +
  geom_hline(yintercept = 1,color="black",alpha=0.3,linetype="dashed") +
  #
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=0), position=position_dodge(.7)) +
  geom_point(position=position_dodge(.7), stat="identity",size=2) +
  #
  ylab(bquote(''*h^2*'')) + xlab("") + 
  scale_color_manual(values = met.brewer("Greek",2)) +
  scale_fill_manual(values = met.brewer("Greek",2)) +
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  scale_y_continuous(limits=c(-0.10,1), breaks = seq(0, 1, by = 0.20)) +
  # # scale_x_discrete(labels=h2_brain$region_name2 %>% unique()) +
  scale_x_discrete(limits=rev(levels(h2_brain$region_name2))) +
  facet_grid(. ~adj, scales="free",space="free") +
  coord_flip() +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  NULL

#----------------------------------------------------------------------
# PLOT COGN
# compare each run to main analysis
h2_cogn_comparison <- ggplot(h2_cogn,
                              aes(x=name,y=h2,width=0.7,
                                  color=covs,fill=covs,
                                  shape=run_name
                                  # alpha=sig
                                  )) + 
  #
  geom_hline(yintercept = 0,color="black",alpha=0.3,linetype="dashed") +
  geom_hline(yintercept = 1,color="black",alpha=0.3,linetype="dashed") +
  #
  geom_errorbar(aes(ymin=estimate-1.96*se, ymax=estimate+1.96*se,width=0), position=position_dodge(.7)) +
  geom_point(position=position_dodge(.7), stat="identity",size=2) +
  #
  ylab(bquote(''*h^2*'')) + xlab("") + 
  scale_color_manual(values = met.brewer("Greek",2)) +
  scale_fill_manual(values = met.brewer("Greek",2)) +
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  scale_y_continuous(limits=c(-0.10,1), breaks = seq(0, 1, by = 0.20)) +
  # # scale_x_discrete(labels=h2_brain$region_name2 %>% unique()) +
  scale_x_discrete(limits=rev(levels(h2_cogn$name))) +
  facet_grid(. ~adj, scales="free",space="free") +
  coord_flip() +
  # theme adjustements
  theme_minimal_vgrid(12) +
  panel_border() +
  NULL
#----------------------------------------------------------------------
legend_h2 <- get_legend(  h2_cogn_comparison +
                            # guides(shape = guide_legend(nrow = 1)) +
                            theme(legend.position = "bottom",legend.direction="vertical") )


# combined plot
h2_comparison<-plot_grid(
  h2_cogn_comparison + theme(legend.position="none") + labs(title="Cognitive measures"),
  NULL,
  h2_brain_comparison + theme(legend.position="none") + labs(title="Reading associated brain measures"),
  nrow=1, rel_widths = c(1,0.2,2), labels=c("A","","B") )

h2_plot <-
  plot_grid(h2_comparison,
                    legend_h2,
                    ncol=1, rel_heights=c(1,0.2)
                    )



#----------------------------------------------------------------------
# compare each run to main
# convert to wide format, each run separate
h2_wide<- h2 %>% dplyr::select(-run_name,-estimate) %>%
  pivot_wider(names_from=run,
              values_from =c(subset,n,h2,se,pval,covs,regression,sig) )

h2_wide_main<- h2_wide %>% 
  mutate_all(.,as.character) %>%
  pivot_longer(cols=c(contains("_lm"),contains("_pheno_subset"),contains("ses_covs")), 
               names_to="var_name",values_to="tmp",
               names_repair = "unique") %>%
  mutate(var=strsplit(var_name,"_") %>% sapply("[[",1),
         run=gsub(
           paste0(paste0("^",unique(var),"_"),collapse="|"),
           "",var_name) ) %>%
  dplyr::select(-var_name) %>%
  pivot_wider(names_from=var,values_from=tmp) %>%
  mutate(across(c(
    "n_main","h2_main","se_main","pval_main",
    "n","h2","se","pval"),as.numeric)) # ** error here to convert to numeric**

h2_wide_main <- h2_wide_main %>% 
  mutate(
    run_name=paste0(
      "subset = ",subset,"; regression = ",regression),
    cat=if_else(category!="Cognitive",paste0("Brain (",adj,")"),category) %>% 
      gsub(" \\(-\\)","",.)
  )

# quick and dirty plot
h2main_comp1 <- ggplot(data=h2_wide_main,
       aes(x=h2_main,y=h2,
           color=covs,fill=covs,
           shape=run_name
           )) +
  geom_abline(slope=1,color="black",size=1,alpha=0.5) +
  geom_errorbar(aes(ymin=h2-1.96*se, ymax=h2+1.96*se)) +
  geom_errorbarh(aes(xmin=h2_main-1.96*se_main, xmax=h2_main+1.96*se_main)) +
  geom_point(size=2,alpha=0.5) +
  ylab(bquote(''*h^2*' (comparison run)')) +
  xlab(bquote(''*h^2*'(main analysis)')) +
  scale_color_manual(values = met.brewer("Greek",2)) +
  scale_fill_manual(values = met.brewer("Greek",2)) +
  # scale_y_continuous(limits=c(-0.05,0.6), breaks = seq(0, 1, by = 0.20)) +
  # # scale_x_discrete(labels=h2_brain$region_name2 %>% unique()) +
  # scale_x_discrete(limits=rev(levels(h2_cogn$name))) +
  facet_grid(. ~cat, scales="free",space="free") +
  # theme adjustements
  theme_minimal_grid(12) +
  panel_border() +
  NULL


h2main_comp2<- h2 %>%
  mutate(cat=if_else(category!="Cognitive",paste0("Brain (",adj,")"),category) %>% 
  gsub(" \\(-\\)","",.)
  ) %>% 
  ggplot(data=.,
         aes(x=n,y=se,
             color=covs,fill=covs,shape=run_name
             )) +
 geom_point(size=2,alpha=0.5) +
  ylab(bquote('Standard Error ('*h^2*')')) +
  xlab("N") +
  scale_color_manual(values = met.brewer("Greek",2)) +
  scale_fill_manual(values = met.brewer("Greek",2)) +
  # scale_y_continuous(limits=c(-0.05,0.6), breaks = seq(0, 1, by = 0.20)) +
  # # scale_x_discrete(labels=h2_brain$region_name2 %>% unique()) +
  # scale_x_discrete(limits=rev(levels(h2_cogn$name))) +
  facet_grid(. ~cat, scales="free",space="free") +
  # theme adjustements
  theme_minimal_grid(12) +
  panel_border() +
  NULL


legend_2 <- get_legend(  h2_brain_comparison +
                           # guides(shape = guide_legend(nrow = 1)) +
                           theme(legend.position = "right",
                                 legend.direction="vertical") )

h2main_comp <- plot_grid(
  h2main_comp2 + theme(legend.position="none"),
  h2main_comp1 + theme(legend.position="none"),
  ncol=1,labels=c("C","D"))

h2_plot2<-plot_grid(
  h2main_comp,legend_2,
  nrow=1,rel_widths = c(1,0.3) )

h2_plot <-
  plot_grid(h2_comparison,
            h2_plot2,
            ncol=1, rel_heights=c(0.8,1.2)
  )

# save
ggsave(h2_plot,
       file=paste0(base_dir,"/output/h2_comparison_acrossRuns.png"),
       height=10,width=12)

