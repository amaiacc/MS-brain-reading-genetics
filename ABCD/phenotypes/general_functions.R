# function to make first letter of each word capital:
## from: https://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#---------------------------------------------------
censor =  function(x, fraction=.005, winsor=FALSE){
  if(length(fraction) != 1 || fraction < 0 || fraction > 1){
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction/2, 1-fraction/2), na.rm = T)
  if(winsor==FALSE) {
    x[ x < lim[1] ] <- NA
    x[ x > lim[2] ] <- NA 
  } else {
    x[ x < lim[1] ] <- lim[1]
    x[ x > lim[2] ] <- lim[2]
  }
  return(x)
}

rm_outliers<-function(x,thr=4){
  m<-mean(x,na.rm=TRUE)
  s<-sd(x,na.rm=TRUE)
  w<-which(x>m+thr*s | x<m-thr*s)
  
  x2<-x
  if(length(w)>0){x2[w]<-NA}
  return(x2)
}

summary_vars <- function(d,n_vars=c(),c_vars=c()) {
  summary_n_vars<- apply(d[,n_vars],2,function(x){
    data.frame(
      n=length(!is.na(x)),
      range=paste(min(x)%>% round(.,digits=3),max(x)%>% round(.,digits=3),sep=";") ,
      median=median(x,na.rm=TRUE) %>% round(.,digits=3),
      mean=mean(x,na.rm=TRUE) %>% round(.,digits=3),
      sd=sd(x,na.rm=TRUE) %>% round(.,digits=3),
      skewness=skewness(x,na.rm=T) %>% round(.,digits=3),
      kurtosis=kurtosis(x,na.rm=T)%>% round(.,digits=3)
    )
  }) %>% do.call("rbind",.)  %>% as.data.frame(stringsAsFactors=FALSE)  %>% mutate(variable=rownames(.))
  
  counts_c_vars<-apply(d[,c_vars],2,function(x){
    t<-table(x)
    t2<-paste(names(t),t,sep=" = ",collapse="\n")
    data.frame(
      counts=t2,
      levels=length(unique(x)),
      n=length(!is.na(x))
    )
  }) %>% do.call("rbind",.) %>% as.data.frame(stringsAsFactors=FALSE) %>% mutate(variable=rownames(.))
  
  summary_vars<-merge(summary_n_vars,counts_c_vars,all=TRUE,stringsAsFactors=FALSE,by=c("variable","n"))
  summary_vars<-summary_vars %>% select(variable,n,levels,counts,range,median,mean,sd,skewness,kurtosis) %>% arrange(levels,variable)
  
  summary_vars[is.na(summary_vars)]<-"-"
  
  return(summary_vars)
}
#---------------------------------------------------