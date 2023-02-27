# load required packages
library(mediation)
library(lme4)
library(dplyr); library(tidyr)

# # trying to edit number of digits for mediation
# as.list(body(mediation:::print.summary.mediate))[11]
# trace(mediation:::print.summary.mediate,
#       at=11,
#       tracer=quote({
#         printCoefmat<-function(x,digits){
#           p<-x[,4]
#           x[,1:3]<-sprintf("%.6f",x[,1:3])
#           x[,4]<-sprintf("%.3f",x[,4])
#           print(x,quote=FALSE,right=TRUE)
#         }
#       }),
#       print=FALSE)
# mediation:::print.summary.mediate(summary(results))

#
run_mediation<-function(dv, iv, mediator, covs, random, d){
  model_covs<-paste("~",paste(c(random,covs),collapse=" + "))
  
  # make sure class is numeric (not matrix as output from the scaling) otherwise mediation complains
  d[,iv]<-as.numeric(d[,iv])
  d[,mediator]<-as.numeric(d[,mediator])
  #---------------------------------------------------
  ## define models
  
  
  ## 1. Path C: total effect: DV  ~ IV
  #  independent variable (iv) dependent variable (dv): reading  ~  CSAstg
  model1<-paste0(dv,model_covs," + ",iv)
  # model1<-paste0(dv ,"~", iv)
  # model1<-paste0(dv,"~",paste(c(covs,iv),collapse=" + "))
  
  ## 2.  Path A: mediator ~ IV
  # The effect of the IV onto the mediator
  ## To establish any mediation, the independent variable (iv) must significantly affect the mediator.
  ## This makes sense, as, for a mediation to take place, the iv must significantly impact the mediator.
  model2<-paste0(mediator,model_covs," + ",iv)
  # model2<-paste0(mediator,"~",iv)
  # model2<-paste0(mediator,"~",paste(c(covs,iv),collapse=" + "))
  
  ## 3. Path C': DV ~ mediator + iv
  # The effect of the mediator on the dependent variable
  # The third step confirms that the mediator affects the dependent variable while controlling for the independent variable. This means, that for a mediation to take place, the mediator must explain more or other parts of the variance in the dependent variable than the independent variable.
  model3<-paste0(dv,model_covs," + ",iv," + ",mediator)
  # model3<-paste0(dv ,"~", iv ,"+", mediator)
  # model3<-paste0(dv,"~",paste(c(covs,iv,mediator),collapse=" + "))
  #---------------------------------------------------
  ## Run models
  m1<-do.call("lmer",list (as.formula(model1),data=d))
  m2<-do.call("lmer",list (as.formula(model2),data=d))
  m3<-do.call("lmer",list (as.formula(model3),data=d))
  #
  # m1<-do.call("lm",list (as.formula(model1),data=d))
  # m2<-do.call("lm",list (as.formula(model2),data=d))
  # m3<-do.call("lm",list (as.formula(model3),data=d))
  # #---------------------------------------------------
  # check models
  anova(m1)
  anova(m2)
  anova(m3)
  #---------------------------------------------------
  ## 4. Causal mediation analysis: calculate the entire model in one
  results = mediation::mediate(m2, m3, treat=iv, 
                               mediator=mediator,
                               covariates=covs,
                               boot=FALSE, # cannot use boot=TRUE with lmer objects
                               # robustSE=TRUE,
                               sims=10000) 
  # boot=FALSE)
  # t<-psych::mediate(y=dv,x=iv,m=mediator,data=d)
  #---------------------------------------------------
  return(results)
}

#---------------------------------------------------
get_mediation_table<-function(x2){
  x<-summary(x2)
  
  #get dv info
  f<-x2$model.y %>% summary()
  f<-f$call$formula %>% as.character() %>% unlist()
  dv<-f[2]
  rm(f)
  # run info
  iv<-x$treat
  med<-x$mediator
  covs<-x$covariates
  
  # get A and B effects
  ## A is the effect of the IV onto the mediator: mediator ~ IV
  A_model<-x2$model.m %>% summary() %>% coef() %>% data.frame()
  A_model$var<-rownames(A_model)
  A_ci<-x2$model.m %>% confint() %>% data.frame()
  A_ci$var<-rownames(A_ci)
  A<-merge(A_model,A_ci,by="var") %>% mutate(p=NA)
  Ap<-A %>% filter(var==iv) %>% dplyr::select(Estimate,X2.5..,X97.5..,p) 
  colnames(Ap)<-c("Estimate","CIlower","CIupper","p")
  rm(A_model,A_ci,A)
  
  ## B is the effect of the mediator on the DV, adjusted for the IV: DV ~ mediator + IV
  B_model<-x2$model.y %>% summary() %>% coef() %>% data.frame()
  B_model$var<-rownames(B_model)
  B_ci<-x2$model.y %>% confint() %>% data.frame()
  B_ci$var<-rownames(B_ci)
  B<-merge(B_model,B_ci,by="var") %>% mutate(p=NA)
  Bp<-B %>% filter(var==med) %>% dplyr::select(Estimate,X2.5..,X97.5..,p)
  Cdp<-B %>% filter(var==iv) %>% dplyr::select(Estimate,X2.5..,X97.5..,p) # C
  colnames(Bp)<-colnames(Cdp)<-c("Estimate","CIlower","CIupper","p")
  rm(B_model,B_ci,B)
  #-----------------------
  # stats from mediation
  d<-cbind(Estimate=x$d.avg,CIlower=x$d.avg.ci[1],CIupper=x$d.avg.ci[2],p=x$d.avg.p)  # ACME
  z<-cbind(Estimate=x$z.avg,CIlower=x$z.avg.ci[1],CIupper=x$z.avg.ci[2],p=x$z.avg.p) # ADE
  t<-cbind(Estimate=x$tau.coef,CIlower=x$tau.ci[1],CIupper=x$tau.ci[2],p=x$tau.p) # Total
  prop<-cbind(Estimate=x$n.avg,CIlower=x$n.avg.ci[1],CIupper=x$n.avg.ci[2],p=x$n.avg.p) # Proportion
  
  m<-rbind(Ap,Bp,d,Cdp,z,t,prop) %>% as.data.frame()
  m$stat<-c("A","B","ACME","C'","ADE","Total","Proportion")
  
  # order
  m<- m[,c("stat","Estimate","CIlower","CIupper","p")]

  # add run info
  m$dv<-dv
  m$iv<-iv
  m$mediator<-med
  m$covariates<-paste(covs,collapse=",")
  
  # clean rownames
  rownames(m)<-1:NROW(m)
  # clean intermediate objects
  rm(Ap,Bp,d,Cdp,z,t,prop)
  rm(dv,iv,med,covs)
  rm(x)
  #
  return(m)
}
