# functions to plot genetic correlations within atlas
corrplot_rg_atlas<-function(atlas=rg_atlas,lr=rg_lr){
  tmp<-subset(atlas,parcellation==parc & measure==m)
  l<- tmp %>% filter(hemisphere=="L") %>% select(region.p1,region.p2,rg,p,hemisphere) # upper triangle: L
  r<- tmp %>% filter(hemisphere=="R") %>% select(region.p2,region.p1,rg,p,hemisphere) # lower triangle: R
  notlat<- tmp %>% filter(is.na(hemisphere)) %>% select(region.p2,region.p1,rg,p,hemisphere)
  lr<-subset(lr,parcellation==parc& measure==m) %>% mutate(region2=region,hemisphere="LR") %>%
    select(region,region2,rg,p2,hemisphere) %>% rename(p=p2) # get pval for different to 1
  colnames(l)<-colnames(r)<-colnames(lr)<-c("p1","p2","rg","p","hemisphere")
  t<-rbind(l,r,lr);rm(l,r,lr)
  t<-subset(t,!is.na(rg))
  # if(sum(t$rg>1,na.rm=TRUE)>0){
  #   t$rg[t$rg>1]<-1 # keep range of rg between -1 and 1 (otherwise corrplot complains)
  # }
  t_rg<-t %>% select(p1,p2,rg) %>% spread(p1,rg) %>% select(-1) %>% as.matrix()
  t_p<-t %>% select(p1,p2,p) %>% spread(p1,p) %>% select(-1) %>% as.matrix()
  rownames(t_rg)<-rownames(t_p)<-colnames(t_rg)
  
  corrplot(corr=t_rg, p.mat=t_p,
           method="circle",
           order="hclust",
           # style
           col=brewer.pal(n=8, name="PuOr"),
           insig="blank",
           is.corr=FALSE,
           addCoef.col=TRUE,
           addCoefasPercent=TRUE,
           type="full",
           diag=TRUE,
           tl.col="black",
           tl.srt=45,
           sig.level=0.05,
           mar=c(0,0,1,0),
           title=paste(m, " (",parc,")", sep=""))
  # return(c)
}
# for not lateralized measures
corrplot_rg_notlat_atlas<-function(atlas=rg_atlas,lr=rg_lr,parc="FAST",m="VOLUME"){
  tmp<-subset(atlas,parcellation==parc & measure==m)
  t1<- tmp %>% filter(is.na(hemisphere)) %>% select(region.p2,region.p1,rg,p)
  colnames(t1)<-c("p1","p2","rg","p")
  t2<-t1 %>% select(p2,p1,rg,p)
  colnames(t2)<-c("p1","p2","rg","p")
  t<-rbind(t1,t2)
  rm(t1,t2)
  t<-subset(t,!is.na(rg))
  # if(sum(t$rg>1,na.rm=TRUE)>0){
  #   t$rg[t$rg>1]<-1 # keep range of rg between -1 and 1 (otherwise corrplot complains)
  # }
  t_rownames<-t %>% select(p1,p2,rg) %>% spread(p1,rg) %>% select(1) %>% pull()
  t_rg<-t %>% select(p1,p2,rg) %>% spread(p1,rg) %>% select(-1) %>% as.matrix()
  t_p<-t %>% select(p1,p2,p) %>% spread(p1,p) %>% select(-1) %>% as.matrix()
  rownames(t_rg)<-rownames(t_p)<- t_rownames
  
  corrplot(corr=t_rg, p.mat=t_p,
           method="circle",
           # order="hclust",
           # style
           col=brewer.pal(n=8, name="PuOr"),
           insig="blank",
           is.corr=FALSE,
           addCoef.col=TRUE,
           addCoefasPercent=TRUE,
           type="lower",
           diag=FALSE,
           tl.col="black",
           tl.srt=45,
           sig.level=0.05,
           mar=c(0,0,1,0),
           title=paste(m, " (",parc,")", sep=""))
  # return(c)
}
