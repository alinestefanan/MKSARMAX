# y: time series
# S: seasonal periodicity
# tau: quantil, when set 0.5 is the median
# link: "logit", "probit" or "cloglog"

sample.mkarma <- function(n,beta0=c(0.0),phi=c(0.0),theta=c(0.0),PHI=c(0.0),THETA=c(0.0),alpha=10,
                            ar=c(0.0),ma=c(0.0),AR=c(0.0),MA=c(0.0),S=12,
                            tau=0.5,freq=12,link="logit")
{
  n=n
  if(phi[1]!=0 & ar[1]==0)
  {
    ar <- 1:length(phi)
  }
  
  if(theta[1]!=0 & ma[1]==0)
  {
    ma <- 1:length(theta)
  }
  
  if(PHI[1]!=0 & AR[1]==0)
  {
    AR <- 1:length(PHI)
  }
  
  if(THETA[1]!=0 & MA[1]==0)
  {
    MA <- 1:length(THETA)
  }
  
  p <- max(ar)
  q <- max(ma)
  P <- max(AR)
  Q <- max(MA)
  m <- 2*max(p,q,P,Q,S*P,S*Q,S*P+p,S*Q+q,na.rm=T)#tirar 2x?
  X<-matrix(rep(1,(n+m)), nrow=(n+m), ncol=1, byrow=F) #default
    ##funções de ligação
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  #multiplicações no eta
  operator<-function(phi,PHI,ar,AR)
  {
    parameters<-c(phi,PHI)
    index<-c(ar,AR)
    j1<-1
    for(j in ar)
    {
      J1<-1
      for(J in AR)
      {      
        parameters<-c(parameters, -phi[j1]*PHI[J1])
        index<-c(index, (j+J))
        J1<-J1+1
      }
      j1<-j1+1
    }
    z<-c()
    z$parameters<-parameters
    z$index<-index
    return(z)
  }
  
  #gerar amostra modified kumaraswamy reparametrizada
  rmkum_a<-function(n,mu,alpha)
  {
    u<- runif(n,min=.Machine$double.eps)
    if(is.na(mu)) mu<-.Machine$double.eps
    if(mu<.Machine$double.eps) mu<-.Machine$double.eps
    if(mu>0.9999999) mu<-0.9999999#1-.Machine$double.eps
    critical<-(alpha-alpha/mu)
    if(is.na(critical)) critical<- -.Machine$double.eps
    if(is.nan(critical)){critical=709.7827}
    if(critical < -36.04365) critical <- -36.04365
    den.cr=log(1-exp(critical))
    if(is.nan(den.cr)){den.cr=-36.04365}
    beta0_cond<-log(1-tau)/den.cr
    y<- alpha/(alpha-log(1-(1-u)^(1/beta0_cond)))
    for (i in 1:n){
    if(y[i] <= .Machine$double.eps| y[i] >= 0.9999999){
      while(y[i] <= .Machine$double.eps| y[i] >= 0.9999999){
      u<- runif(1,min=.Machine$double.eps)
      y[i]<- alpha/(alpha-log(1-(1-u)^(1/beta0_cond)))
      }
    }
    }
    return(y)
  }
  ar_par_index <- operator(phi,PHI,ar,S*AR)
  ma_par_index <- operator(theta,THETA,ma,S*MA)
  ar_par <- ar_par_index$parameters
  ar_ind <- ar_par_index$index
  
  ma_par <- ma_par_index$parameters
  ma_ind <- ma_par_index$index
 
  ynew <-rep(beta0[1],(n+m)) # primeiro valor aleatório= beta0[1]
  mu <- linkinv(ynew)
  
  error<-rep(0,n+m) # E(error)=0 
  eta<- y <- rep(NA,n+m)
  parte_a<-rep(0,n+m)
  parte_b<-rep(0,n+m)
  parte_c<-rep(0,n+m)
  mu0<-rep(NA,n)
  critical<-c()

  for(i in (m+1):(n+m)) 
  {
    parte_a[i]<- X[i,]%*%as.matrix(beta0)
    parte_b[i]<- sum(ar_par*(ynew[i-ar_ind]))
    parte_c[i]<- sum(ma_par*error[i-ma_ind])
    
    eta[i]  <- parte_a[i] +  parte_b[i] - parte_c[i]

    mu[i]   <- linkinv(eta[i])

    y[i]    <- rmkum_a(1, mu[i], alpha)###gerando amostra modified kumaraswamy reparametrizada
   
    ynew[i] <- linkfun(y[i])

    error[i] <- ynew[i]-eta[i] #residuals on the predictor scale
  }

  return( ts(y[(m+1):(n+m)],frequency=freq) )
}
