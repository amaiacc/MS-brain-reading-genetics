# Edited to take any input, specified as an argument
# From: https://github.com/ENIGMA-git/ENIGMA/tree/master/Genetics/ENIGMA1/Imputation
#-------------------------------------
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
#-------------------------------------
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
# args<-c("ABCD_release_3.0_QCed_pruned_1KG_mds2R.mds","ABCD_release_3.0_QCed_1KG.pruned","ABCD_release_3.0_QCed_pruned_1KG_mds","ABCD_release_3.0_QCed.pop_strat_outliers.outliers","ABCD_release_3.0_QCed")

print(args)
#-------------------------------------
# get population codes: https://www.internationalgenome.org/category/population/
popcodes<-read.csv("PopulationCodes.csv",stringsAsFactors=FALSE)
#'s/ASW/3/g' -e 's/CEU/4/g' -e 's/CHB/5/g' -e 's/CHS/6/g' -e 's/CLM/7/g' -e 's/FIN/8/g' -e 's/GBR/10/g' -e 's/IBS/11/g' -e 's/JPT/12/g'
#' -e 's/LWK/13/g' -e 's/MXL/14/g' -e 's/PUR/15/g' -e 's/TSI/16/g' -e 's/YRI/17/g'
# create population correspondences
pops<-rbind(c("ASW","3"),
            c("CEU","4"),
            c("CHB","5"),
            c("CHS","6"),
            c("CLM","7"),
            c("FIN","8"),
            c("GBR","10"),
            c("IBS","11"),
            c("JPT","12"),
            c("LWK","13"),
            c("MXL","14"),
            c("PUR","15"),
            c("TSI","16"),
            c("YRI","17"),
            c("Dataset","-9"))
pops<-data.frame(pops)
colnames(pops)<-c("Population.Code","Population")
#
pops<-merge(pops,popcodes,all.x=TRUE,stringsAsFactors=FALSE)
pops$Super.Population.Code<-as.character(pops$Super.Population.Code)
pops<-pops[order(pops$Super.Population.Code),]
pops$Population.Code[is.na(pops$Population.Code)]<-"Dataset"
#-------------------------------------
# read data

## MDS
mds.cluster <- read.table(args[1],header=TRUE,sep=" ",stringsAsFactors=FALSE) #read in file as a dataframe

# Read corresponding .fam file, to get population
infile<-args[2] # specify input file for .fam
fam<-read.table(paste(infile,".fam",sep=""),header=FALSE,sep=" ")
colnames(fam)<-c("FID","IID","pID","mID","sex","Population")

# Specify output file
outfile<- args[3] 

# PopStr outliers, if file exists
if (file.exists(pattern=paste(args[4],sep=""))){

  outliers<-read.table(args[4])
  colnames(outliers)<-c("FID","IID")
  outliers$outlier<-1 # flag outliers
  
  fam<-merge(fam,outliers,all=TRUE)
  if (length(is.na(fam$outlier))>1) { fam$outlier[is.na(fam$outlier)]<-0 }
  fam$outlier<-as.factor(fam$outlier)

} else {
  fam$outlier<-NA
}

# add info about batches
r<-args[5]
if (file.exists(pattern=paste(r,".batches",sep=""))){
  batch<-read.table(paste(r,".batches",sep=""),header=TRUE)
  fam<-merge(fam,batch,by=c("FID","IID"),all.x=TRUE)
  rm(batch)
} else {
  fam$batch<-NA
  fam$batch[which(fam$Population==100|fam$Population==-9)]<-args[5]
}
rm(r)
#-------------------------------------
# Merge with mds file
mds.cluster.fam<-merge(mds.cluster,fam,all=TRUE)
mds.cluster.fam2<-merge(mds.cluster.fam,pops,all.x=TRUE)

mds.cluster.fam2$Population.Code[which(!is.na(mds.cluster.fam2$batch))]<-mds.cluster.fam2$batch[which(!is.na(mds.cluster.fam2$batch))]

# recode levels
# lvls<-unique(as.character(pops$Population.Code) )
lvls<-c(as.character(pops$Population.Code)[-length(unique(pops$Population.Code))],unique(fam$batch)[!is.na(unique(fam$batch))])
mds.cluster.fam2$Population.Code <- factor(mds.cluster.fam2$Population.Code, levels = lvls)
# define colours for ethnicity, given Super Population Code
mycols<-c(brewer.pal(4, "Reds")[-1],brewer.pal(4, "Purples")[-1],brewer.pal(4, "Greens")[-1],brewer.pal(6, "Blues")[-1])
mycols<-c(mycols,brewer.pal((length(lvls)-length(mycols))+1,"Greys")[-1])
rm(lvls)

#-------------------------------------
# Plot Dimensions 1 to 3
#-------------------------------------
for (x in 1:3){
  pcx=paste0("C",x)
  pcy=paste0("C",x+1)
  p_1kg <- ggplot(data=subset(mds.cluster.fam2,(is.na(batch))), aes_string(x=pcx,y=pcy,color="Population.Code")) + 
    geom_point(alpha=0.3,size=2) +
    scale_colour_manual(values = c(mycols) ) +
    NULL
  
  p_1kg_d <- p_1kg + 
              geom_point(data=subset(mds.cluster.fam2,!(is.na(batch))),aes(shape=batch),color="black",alpha=0.3) + # ,fill=outlier,size=outlier
              # Fix legends to be informative
              scale_shape_manual(name = "Dataset", values = c(21,24)) +
              scale_fill_manual(name = "Outlier", values = c("white","black")) +
              # defining size with 2 marginally different values
              scale_size_manual(name = "Population stratification", values = c(2, 2.01),labels=c("","Outlier")) +
              # Remove fill legend and replace the fill legend using the newly created size
              guides(fill = "none", 
              size = guide_legend(override.aes = list(shape = c(1, 16)))) +
              theme(legend.position = "bottom") +
              NULL
  # zoom in to the sample cluster region (without outliers)
  if (sum(!is.na(mds.cluster.fam2$outlier)>0)){
    x_range<-range(subset(mds.cluster.fam2,!(is.na(batch))&outlier==0)[,pcx])
    y_range<-range(subset(mds.cluster.fam2,!(is.na(batch))&outlier==0)[,pcy])
    p_1kg_d_zoom <- p_1kg_d + 
                      coord_cartesian(xlim = x_range,ylim= y_range) +
                      # adjust legends to remove non-informative info
                      guides(fill="none",size="none") +
                      theme(legend.position="bottom") + 
                      NULL
                      
  }  
  # cairo_pdf(
  # pdf(file=paste(outfile,"plot_",pcx,"vs",pcy,".pdf",sep=""),width=12) #Use this file to look at your MDS plot if you do not
  png(file=paste(outfile,"plot_",pcx,"vs",pcy,".png",sep=""),width=800,height=800) #Use this file to look at your MDS plot if you do not
  print(p_1kg_d)
  dev.off()
  
  
  if (exists("p_1kg_d_zoom")){
    # cairo_pdf(
    # pdf(file=paste(outfile,"plot_",pcx,"vs",pcy,"_zoomed.pdf",sep=""),width=12) #Use this file to look at your MDS plot if you do not 
    png(file=paste(outfile,"plot_",pcx,"vs",pcy,"_zoomed.png",sep=""),width=800,height=800) #Use this file to look at your MDS plot if you do not
    print(p_1kg_d_zoom)
    dev.off()
  }
}



#-------------------------------------
