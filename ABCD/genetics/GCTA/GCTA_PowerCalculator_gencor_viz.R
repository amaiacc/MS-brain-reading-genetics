#----------------------------------------------------------------------
# GCTA power calculator - for genetic correlation
## between a behavioural trait/ brain measure
## given h2 estimates and N
#----------------------------------------------------------------------
library(tidyr)
library(dplyr)

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
#----------------------------------------------------------------------

#----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
# define pattern for this run
root="ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned_rm05_adj"

# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
working_dir=paste(dir,"resources/datasets/ABCD/data/working_data/genotyping/GCTA/",sep="")
scripts_dir=paste(dir,"resources/datasets/ABCD/scripts/genetics/GCTA/",sep="")
out_dir=paste(working_dir,"/output/",sep="")
setwd(out_dir)
#
source(paste(scripts_dir,"../../phenotypes/color_functions.R",sep=""))
#-----------------------------------------------
resultsQT<-read.csv(file=paste0(out_dir,"power_calc_results_",root,".csv",sep=""))
resultsQT$hemisphere[which(is.na(resultsQT$hemisphere))]<-"NotLat"
resultsQT$hemisphere<-factor(resultsQT$hemisphere,levels=c("L","R","NotLat"))
resultsQT$category=resultsQT$category.1
#
# clean region name
resultsQT$region<-gsub("_adj|global|Global","",resultsQT$region)
resultsQT$region_name <- resultsQT$region %>% gsub("smri_thick_cort.|smri_area_cort.|destrieux_","",.) %>% 
  gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% gsub("AREA_|area_|THICK_|thick_","",.) %>% 
  toupper()
resultsQT$measure<-NA
resultsQT$measure[grep("area",resultsQT$region)]<-"AREA"
resultsQT$measure[grep("thick",resultsQT$region)]<-"THICKNESS"
# 
resultsQT<-subset(resultsQT,hemisphere=="L")

#-----------------------------------------------
# Make plots
#-----------------------------------------------
# colors_t<-brewer.pal(length(unique(resultsQT$region[which(resultsQT$measure=="THICKNESS")]))+1, "Set1")[-6]
# colors_a<-brewer.pal(length(unique(resultsQT$region[which(resultsQT$measure=="AREA")])), "Set2")
# colors<-c(colors_a,colors_t)
braincolor_codes<-colors_brain(data=resultsQT,atlas="destrieux",p_col="trait2")

resultsQT<-merge(resultsQT,braincolor_codes,by.x=c("region","measure","hemisphere"),
                 by.y=c("Region","measure","hemisphere")) %>% arrange(region)
resultsQT$region<-factor(resultsQT$region,
                       levels=levels(braincolor_codes$Region)) 
resultsQT<-resultsQT %>% arrange(region) %>% mutate(region_name=factor(region_name,levels=unique(region_name)))

resultsQT$measure[grep("area",resultsQT$region)]<-"CSA"
resultsQT$measure[grep("thick",resultsQT$region)]<-"CT"


colors_a <- resultsQT %>% filter(measure=="AREA"&hemisphere=="L") %>% pull(color) %>% unique()
colors_t <- resultsQT %>% filter(hemisphere=="L" & measure=="THICKNESS") %>% pull(color)  %>% unique()
colors <- resultsQT  %>% filter(!is.na(hemisphere)) %>% pull(color) %>% unique()

plot_power<-function(x,colors){ 
  x2<-gsub("nihtbx_|_uncorrected","",x)
  d<-subset(resultsQT,trait1==x)
  

  gl <- ggplot(data=d %>% filter(category=="Brain global"),aes(x=rg,y=pwr,color=region_name)) +
    geom_hline(yintercept =0.8,color="black",linetype="dashed") +
    geom_line(size=2) +
    # scale_color_brewer(palette="Set2") +
    scale_color_manual(values=colors) +
    ylim(c(0,1)) +
    labs(x = bquote(''*rho*'(G)'), y = "Power", title = "") +
    theme_minimal_grid(12) +
    panel_border() +
    guides(colour=guide_legend(nrow=3)) +
    theme(legend.position = "bottom") +
    NULL
  
  rs <- ggplot(data=d,aes(x=rg,y=pwr,color=region_name)) + 
    geom_hline(yintercept =0.8,color="black",linetype="dashed") + 
    geom_point(size=1.5,alpha=0.5) +
    geom_line(size=1) +
    scale_color_manual(values=colors,name="") +
    ylim(c(0,1)) +
    facet_grid(measure ~ adj ,scales = "free_x") +
    labs(x = bquote(''*rho*'(G)'), y = "Power", title = "") + 
    theme_minimal_grid(12) +
    panel_border() +
    guides(colour=guide_legend(nrow=3)) +
    theme(legend.position = "bottom") +
    NULL
  
  # plot_row<-plot_grid(gl + theme(legend.position = "bottom"),
  #                     rs + theme(legend.position = "bottom"),rel_widths = c(1,2))
  
  plot_row<-rs
  
  legend0 <- get_legend(gl + scale_color_brewer(palette="Set1",guide="none"))
  legend1 <- get_legend(gl + scale_linetype(guide="none") )
  legend2 <- get_legend(rs + scale_linetype(guide="none") )
  
  bottom_row <- plot_grid(
    legend1, legend0, legend2,
    nrow = 2
  )
  
  # now add the title
  title <- ggdraw() + draw_label(
      paste0("Power to detect genetic correlations with ",x2),
      fontface = 'bold',
      x = 0,hjust = 0 ) +
    theme( # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)  )
  
  m="GREML"
  subtitle <- ggdraw() + draw_label(
    paste0("Method: ",m),
    x = 0,hjust = 0 ) +
    theme( # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)  )
  p<-plot_grid(
    title, subtitle, 
    plot_row,
    # bottom_row,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05,0.05, 1)
  )
  return(p)
}

#-----------------------------------------------
# plot 
mytheme<-  theme(axis.title.x=element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1,face="bold"),
                 axis.text=element_text(size=12),
                 legend.position="bottom",legend.direction = "vertical")
p<-plot_power(x = "nihtbx_reading_uncorrected",colors=colors) + guides(colour=guide_legend(nrow=3))
ggsave(p,file=paste0(out_dir,"power_rg_reading_ROIs.png"),width=8,height=7)

# plot_power(x = "nihtbx_picvocab_uncorrected",colors=colors)
# plot_power(x = "nihtbx_cryst_uncorrected",colors=colors)
# plot_power(x = "nihtbx_fluidcomp_uncorrected",colors=colors)
#-----------------------------------------------
