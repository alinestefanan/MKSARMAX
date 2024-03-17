# y: time series
# S: seasonal periodicity
# exvar: covariate column matrix
# tau: quantil, when set 0.5 is the median
# resid: 1 = quantile residual, 2 = deviance residual
# link: "logit", "probit" or "cloglog"
# steps: how many steps to forecast

mkarma <- function(y,ar=c(0.0),ma=c(0.0),AR=c(0.0),MA=c(0.0),S=12,exvar=matrix(NA, nrow=1, ncol=1, byrow=F),
                tau=0.5,resid=1,aclag=10,steps=12,validation=T,graph=T,print=F,check=F,link="logit")
  {
  if (ar[1]!=0 && ma[1]!=0 && AR[1]!=0 && MA[1]!=0){
    if(#is.matrix(exvar) && is.na(exvar[1,1])==F ||
       is.na(exvar[1])==F)
    {
      source("mkarmaSARMAX.R")
      EMV.mkarmaSARMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkarmaSARMA.R")
    EMV.mkarmaSARMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
  }#armaSARMA
  else if (ar[1]!=0 && ma[1]!=0 && AR[1]==0 && MA[1]==0){
    if(#is.matrix(exvar) && is.na(exvar[1,1])==F || 
       is.na(exvar[1])==F)
    {
      source("mkarmaX.R")
      EMV.mkarmaX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
      source("mkarma.R")
      EMV.mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
  }#arma
  else if (ar[1]!=0 && ma[1]==0 && AR[1]==0 && MA[1]==0){ 
    if(#is.matrix(exvar) && is.na(exvar[1,1])==F || 
       is.na(exvar[1])==F)
    {
      source("mkarX.R")
      EMV.mkarX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
      source("mkar.R")
      EMV.mkar(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
  }#ar
  else if (ar[1]==0 && ma[1]!=0 && AR[1]==0 && MA[1]==0){ 
    if(#is.matrix(exvar) && is.na(exvar[1,1])==F || 
       is.na(exvar[1])==F)
    {
      source("mkmaX.R")
      EMV.mkmaX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
      source("mkma.R")
      EMV.mkma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
  }#ma
  else if (ar[1]!=0 && ma[1]==0 && AR[1]!=0 && MA[1]==0){ 
    if( is.na(exvar[1])==F)
    {
      source("mkarSARX.R")
      EMV.mkarSARX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkarSAR.R")
    EMV.mkarSAR(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
     }#arSAR
  else if (ar[1]==0 && ma[1]==0 && AR[1]!=0 && MA[1]!=0){
    if( is.na(exvar[1])==F)
    {
      source("mkSARMAX.R")
      EMV.mkSARMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkSARMA.R")
    EMV.mkSARMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
     }#SARMA
  else if (ar[1]==0 && ma[1]==0 && AR[1]!=0 && MA[1]==0){
    if( is.na(exvar[1])==F)
    {
      source("mkSARX.R")
      EMV.mkSARX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkSAR.R")
    EMV.mkSAR(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#SAR
  else if (ar[1]!=0 && ma[1]==0 && AR[1]!=0 && MA[1]!=0){
    if( is.na(exvar[1])==F)
    {
      source("mkarSARMAX.R")
      EMV.mkarSARMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkarSARMA.R")
    EMV.mkarSARMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#arSARMA
  else if (ar[1]!=0 && ma[1]!=0 && AR[1]==0 && MA[1]!=0){
    if( is.na(exvar[1])==F)
    {
      source("mkarmaSMAX.R")
      EMV.mkarmaSMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkarmaSMA.R")
    EMV.mkarmaSMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#armaSMA
  else if (ar[1]==0 && ma[1]==0 && AR[1]==0 && MA[1]!=0){
    if( is.na(exvar[1])==F)
    {
      source("mkSMAX.R")
      EMV.mkSMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkSMA.R")
    EMV.mkSMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#SMA
  else if (ar[1]!=0 && ma[1]!=0 && AR[1]!=0 && MA[1]==0){
    if( is.na(exvar[1])==F)
    {
      source("mkarmaSARX.R")
      EMV.mkarmaSARX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkarmaSAR.R")
    EMV.mkarmaSAR(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#armaSAR
  else if (ar[1]==0 && ma[1]!=0 && AR[1]!=0 && MA[1]!=0){
    if( is.na(exvar[1])==F)
    {
      source("mkmaSARMAX.R")
      EMV.mkmaSARMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkmaSARMA.R")
    EMV.mkmaSARMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#maSARMA
  else if (ar[1]!=0 && ma[1]==0 && AR[1]==0 && MA[1]!=0){
    if( is.na(exvar[1])==F)
    {
      source("mkarSMAX.R")
      EMV.mkarSMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkarSMA.R")
    EMV.mkarSMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#arSMA
  else if (ar[1]==0 && ma[1]!=0 && AR[1]!=0 && MA[1]==0){
    if( is.na(exvar[1])==F)
    {
      source("mkmaSARX.R")
      EMV.mkmaSARX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkmaSAR.R")
    EMV.mkmaSAR(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#maSAR
  else if (ar[1]==0 && ma[1]!=0 && AR[1]==0 && MA[1]!=0){
    if( is.na(exvar[1])==F)
    {
      source("mkmaSMAX.R")
      EMV.mkmaSMAX(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,exvar=exvar,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }else{
    source("mkmaSMA.R")
    EMV.mkmaSMA(y=y,ar=ar,ma=ma,AR=AR,MA=MA,S=S,tau=tau,resid=resid,aclag=aclag,steps=steps,validation=validation,graph=graph, print=print, check=check,link=link)
    }
    }#maSMA
  else {
    warning("NO ARMA STRUCTURE INFORMED, MODIFIED KUMARASWAMY DISTRIBUTION ADJUSTED CONSIDERING CONSTANT MEDIAN")
    #função densidade de probabilidade modified kumaraswamy
    dmk <- Vectorize(function(y,alpha,beta,log = FALSE){
      logden <- log((alpha*beta*exp(alpha-alpha/y)*(1-exp(alpha-alpha/y))^(beta-1))/(y^2))
      val <- ifelse(log, logden, exp(logden)) 
      return(val)
    }) #função Murilo
    dmk_alpha<-function(alpha){
      beta=log(1-tau)/log(1-exp(alpha-alpha/median(y)))
      density<-(alpha*beta*exp(alpha-alpha/y)*(1-exp(alpha-alpha/y))^(beta-1))   /(y^2)
      density[is.na(density)] <- .Machine$double.eps
      density[is.nan(density)] <- .Machine$double.eps
      #if(is.na(density)) density<-.Machine$double.eps
      #if(is.nan(density)) density<-.Machine$double.eps
      for (i in 1:length(density)){
        if(density[i]<.Machine$double.eps) density[i]<-.Machine$double.eps
      }
      return(density)
    }
    loglik<-function(z)
    {
      mu<- z[1]
      alpha<-z[2]
      ################################################################################################
      ll <- suppressWarnings(sum(dmk(y, alpha, log(1-tau)/log(1-exp(alpha-alpha/mu)), log = TRUE))) #log-density modified kumaraswamy quantile re-parametrization
      ################################################################################################
      sum(ll)
    }
    escore<-function(z)
    {
      mu<- z[1]
      alpha<-z[2]
    ####INÍCIO DERIVADA DA LOGVEROSSIMILHANÇA EM RELAÇÃO A MU
    ###########################################################################################################
    num1<-exp(alpha)*alpha*(log(1-exp(alpha-alpha/mu))+log(1-tau)*log(1-exp(alpha-alpha/y)))
    den1<-(mu^2)*(exp(alpha/mu)-exp(alpha))*((log(1-exp(alpha-alpha/mu)))^2)
    mustar<-sum(num1/den1)
    ########################################################################################################### 
    ####FIM DERIVADA DA LOGVEROSSIMILHANÇA EM RELAÇÃO A MU
    
    ##############################################################################
    ##### INÍCIO DERIVADA DA LOGVEROSSIMILHANÇA EM RELAÇÃO A alpha
    num1<-(mu-1)*exp(alpha)*log(1-tau)*log(1-exp(alpha-alpha/y))
    den1<- mu*(exp(alpha/mu)-exp(alpha))*(log(1-exp(alpha-alpha/mu)))^2
    num2.part1<-(y-1)*exp(alpha-alpha/y)
    num2.part2<-(log(1-tau)/log(1-exp(alpha-alpha/mu))-1)
    den2<-y-y*exp(alpha-alpha/y)
    den3<- mu*(exp(alpha/mu)-exp(alpha))*log(1-exp(alpha-alpha/mu))
    Ualpha<-sum(1+1/alpha+num1/den1-num2.part1*num2.part2/den2+((mu-1)*exp(alpha))/den3-1/y)
    ##### FIM DERIVADA DA LOGVEROSSIMILHANÇA EM RELAÇÃO A alpha
    ##############################################################################
    return(c(mustar,Ualpha))
    }
    mu<-median(y)
    library(GenSA)
    on.dmk.alpha<-function(alpha){-sum(log(dmk_alpha(alpha)))}
    
    gen.semchute <- GenSA(lower = c(.Machine$double.eps), 
                          upper = c(100),
                          fn = on.dmk.alpha, control=list(max.time=2))
    alpha<-gen.semchute$par
    par_mu_alpha<-c(mu,alpha)
    opt<-optim(par_mu_alpha, fn=loglik, gr=escore,method = "BFGS", control = list(fnscale=-1))
    #plot(y,type="l")
    #plot(dmk(y, opt$par[2], log(1-tau)/log(1-exp(opt$par[2]-opt$par[2]/opt$par[1])), log = TRUE),type = "l")
    hist(y,prob=T)
    x<-seq(0,1,0.01)
    curve(dmk(x, opt$par[2], log(1-tau)/log(1-exp(opt$par[2]-opt$par[2]/opt$par[1])), log = F),add=T)
    result1 <- matrix(round(c(opt$par),4), nrow=length(opt$par), ncol=1, byrow=F)
    colnames(result1) <- c("EMV")
    rownames(result1)<-c("median", "alpha estimator")
    print(result1)
    result2<-matrix(round(c(opt$value),4),nrow=1,ncol=1)
    rownames(result2)<-c("")
    colnames(result2) <-c("Log-likelihood") 
    print(result2)
    }
}



