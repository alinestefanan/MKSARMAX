# y: time series
# S: seasonal periodicity
# exvar: covariate column matrix
# tau: quantil, when set 0.5 is the median
# resid: 1 = quantile residual, 2 = deviance residual
# link: "logit", "probit" or "cloglog"
# steps: how many steps to forecast

EMV.mkarmaX <- function(y,ar=c(0.0),ma=c(0.0),AR=c(0.0),MA=c(0.0),S=12,exvar=matrix(NA, nrow=1, ncol=1, byrow=F),
                        tau=0.5,resid=1,aclag=10,steps=12,validation=T,graph=T,print=T,check=F,link="logit")
{
  if(validation==T){
    n <- length(y)-steps 
  }else{
    n <- length(y) 
  }
  exvar=as.matrix(exvar)
  k=if(is.matrix(exvar)){ncol(exvar)}else{1}
  
  X <- matrix(c(rep(1,(n+steps)),exvar), nrow=(n+steps), ncol=(k+1), byrow=F)
  ##funções de ligação
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == "logit"))
  {  
    stats <- make.link(linktemp)
    dif2link = function(muhat) eval(D(expression(-1/((muhat-1)*muhat)),"muhat"),list(muhat=muhat))
  }else if (any(linktemp == "probit"))
  {  
    stats <- make.link(linktemp)
    library(VGAM)
    dif2link = function(muhat) probitlink(muhat,  deriv = 2)#segunda derivada
  }else if (any(linktemp == "cloglog"))
  {  
    stats <- make.link(linktemp)
    dif2link = function(muhat) eval(D(expression(1/(log(1-muhat)*(muhat-1))),"muhat"),list(muhat=muhat))#segunda derivada
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv,
                         mu.eta = stats$mu.eta,#derivada de mu em relação a eta
                         diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  #função densidade de probabilidade modified kumaraswamy
  dmk <- Vectorize(function(y,alpha,beta,log = FALSE){
    critical.y=exp(alpha-alpha/y)
    for (i in 1:length(y)){
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
    }
    
    density<-(alpha*beta*critical.y*(1-critical.y)^(beta-1))   /(y^2)
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    for (i in 1:length(density)){
      if(density[i]<.Machine$double.eps) density[i]<-.Machine$double.eps
      if(density[i]>.Machine$double.xmax) density[i]<-.Machine$double.xmax
    }
    logden <- log(density)
    val <- ifelse(log, logden, exp(logden)) 
    return(val)
  }) #função Murilo
  
  #modified kumaraswamy cumulative distribution function 
  pmk <- Vectorize(function(q,alpha,beta,log.p = FALSE){
    cdf=rep(0,length(q))
    for(i in length(q)){
      cdf[i] <-  1-(1-exp(alpha-alpha/q[i]))^beta[i]
      if (is.na(cdf[i])) cdf[i]<-.Machine$double.eps
      if (cdf[i]<.Machine$double.eps) cdf[i]<-.Machine$double.eps
      if (cdf[i]>0.9999999) cdf[i]<-0.9999999#1-.Machine$double.eps
    }
    val <- ifelse(log.p, log(cdf), cdf)
    return(val)
  })#função Murilo
  
  dmk_alpha<-function(alpha){
    beta=log(1-tau)/log(1-exp(alpha-alpha/median(y1)))
    if (is.na(beta)){beta=.Machine$double.eps}
    if(beta<.Machine$double.eps) beta<-.Machine$double.eps
    if (is.infinite(beta)){beta=.Machine$double.xmax}
    critical.y=exp(alpha-alpha/y1)
    for (i in 1:length(y1)){
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
    }
    density<-(alpha*beta*critical.y*(1-critical.y)^(beta-1))   /(y1^2)
    density[is.na(density)] <- .Machine$double.eps
    density[is.nan(density)] <- .Machine$double.eps
    for (i in 1:length(density)){
      if(density[i]<.Machine$double.eps) density[i]<-.Machine$double.eps
    }
    return(density)
  }
  
  Monti.test<-function (x, lag = 1, type = c("Ljung-Box"), fitdf = 0) ##Ljung-Box test as the Box.test function using pacf replacing acf
  {
    if (NCOL(x) > 1) 
      stop("x is not a vector or univariate time series")
    DNAME <- deparse(substitute(x))
    type <- match.arg(type)
    cor <- pacf(x, lag.max = lag, plot = FALSE, na.action = na.pass)
    n <- sum(!is.na(x))
    PARAMETER <- c(df = lag - fitdf)
    obs <- cor$acf[1:(lag )]
    if (type == "Ljung-Box") {
      METHOD <- "Monti test via Ljung-Box test"
      STATISTIC <- n * (n + 2) * sum(1/seq.int(n - 1, n - lag) * 
                                       obs^2)
      PVAL <- 1 - pchisq(STATISTIC, lag - fitdf)
    }
    names(STATISTIC) <- "X-squared"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
                   p.value = PVAL, method = METHOD, data.name = DNAME), 
              class = "htest")
  }
  vector.root=function(order,coeff){
    if(sum(order)!=0){
      o=seq(1:max(order))
      p0=rep(0,max(order))
      for (i in 1:length(order)){
        p0[o==order[i]]=coeff[i]
      }
      return(c(1,p0))
    }else{return(c(1,0))}
  }
  p <- max(ar)
  q <- max(ma)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  ynew <-linkfun(y[1:n])     
  
  loglik <- function(z) 
  {
    beta0 <- z[1]
    beta <- z[2:(k+1)]
    phi = z[(k+2):(p1+k+1)] 
    theta =z[(p1+k+2):(p1+q1+k+1)]
    alpha <- z[(p1+q1+k+2)] # alpha parameter
    
    error<-rep(0,n) # E(error)=0 
    eta<-rep(NA,n)
    mu0<-rep(NA,n)
    
    critical<-c()
    den.cr<-c()
    for(i in (m+1):n)
    {
      eta[i] <- X[i,1]*beta0 + X[i,2:ncol(X)]%*%as.matrix(beta) + sum(phi*(ynew[i-ar]-X[i-ar,2:ncol(X)]%*%as.matrix(beta)  ) ) - sum(theta*error[i-ma])
      error[i] <- ynew[i]-eta[i] #residuals on the predictor scale
    }
    
    mu <- linkinv(eta[(m+1):n])
    y1<-y[(m+1):n]
    
    critical<-c()
    den.cr<-c()
    for(i in 1:(n-m))
    {
      if(is.na(mu[i])) mu[i]<-.Machine$double.eps
      if(mu[i]<.Machine$double.eps) mu[i]<-.Machine$double.eps
      if(mu[i]>(0.9999999)) mu[i]<-0.9999999#1-.Machine$double.eps
      
      critical[i]<-(alpha-alpha/mu[i])
      if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
      if(is.nan(critical[i])){critical[i]=-36.04365}
      if(critical[i] < -36.04365) critical[i] <- -36.04365
      
      den.cr[i]=log(1-exp(critical[i]))
      if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
    } 
    l=critical.y=critical.ly=c()
    for (i in 1:length(y1)){
      critical.y[i]=exp(alpha-alpha/y1[i])
      critical.ly[i]<-log(1-critical.y[i])
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
      
      l[i]=log(alpha)+alpha-alpha/y1[i]+(log(1-tau)/den.cr[i] -1)*critical.ly[i]-2*log(y1[i])+log(log(1-tau)/den.cr[i])
    }
    # print("l");print(sum(l))
    #ll <- dmk(y1, alpha, log(1-tau)/den.cr, log = TRUE)#log-density modified kumaraswamy quantile re-parametrization
    #print("ll");print(sum(ll))
    return(sum(l))
  }#fim loglik
  
  score <- function(z) 
  {
    beta0 <- z[1]
    beta <- z[2:(k+1)]
    phi = z[(k+2):(p1+k+1)] 
    theta =z[(p1+k+2):(p1+q1+k+1)]
    alpha <- z[(p1+q1+k+2)] # alpha parameter
    
    error<-rep(0,n) # E(error)=0 
    eta<-rep(NA,n)
    mu0<-rep(NA,n)
    
    critical<-c()
    den.cr<-c()
    for(i in (m+1):n)
    {
      eta[i] <- X[i,1]*beta0 + X[i,2:ncol(X)]%*%as.matrix(beta) + sum(phi*(ynew[i-ar]-X[i-ar,2:ncol(X)]%*%as.matrix(beta)) ) - sum(theta*error[i-ma])
      error[i] <- ynew[i]-eta[i] #residuals on the predictor scale
    }
    
    mu <- linkinv(eta[(m+1):n])
    mu.eta <-  link$mu.eta
    
    critical<-c()
    den.cr<-c()
    for(i in 1:(n-m))
    {
      if(is.na(mu[i]))mu[i]<-.Machine$double.eps
      if(mu[i]<.Machine$double.eps) mu[i]<-.Machine$double.eps
      if(mu[i]>(0.9999999)) mu[i]<-0.9999999#1-.Machine$double.eps
      critical[i]<-(alpha-alpha/mu[i])
      if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
      if(is.nan(critical[i])){critical[i]=-36.04365}
      if(critical[i] < -36.04365) critical[i] <- -36.04365
      den.cr[i]=log(1-exp(critical[i]))
      if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
    }
    y1<-y[(m+1):n]
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    critical.y=exp(alpha-alpha/y1)
    critical.ly<-log(1-critical.y)
    for (i in 1:length(y1)){
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
    }
    num1<-exp(alpha)*alpha*(den.cr+log(1-tau)*critical.ly)
    den1<-(mu^2)*(exp(alpha/mu)-exp(alpha))*((den.cr)^2)
    mustar<-num1/den1
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU      
    mT <- diag(mu.eta(eta[(m+1):n]))
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)# theta
    for(i in 1:(n-m))
    {
      for(j in 1:q1)
      {
        R[i,j] <- -sum(error[i+m-ma[j]])
      }
    }
    
    A <- matrix(rep(NA,(n-m)*(p1)),ncol=(p1))# ar
    for(i in 1:(n-m))
    {
      for(j in 1:p1)
      {
        A[i,j] <- (ynew[i+m-ar[j]]-X[i+m-ar[j],2:ncol(X)]%*%as.matrix(beta))   
      }
    }    
    
    B0 <- matrix(rep(NA,(n-m)),ncol=1)#intercept
    for(i in 1:(n-m))
    {
      B0[i,1] <- X[i+m,1] 
    }
    
    B <- matrix(rep(NA,(n-m)*k),ncol=k)#covariates
    for(i in 1:(n-m))
    {
      for(j in 1:k)
      {
        B[i,j] <- X[i+m,1+j] - sum(phi*X[i+m-ar,1+j])
      }
    }
    deta.dbeta0 <- matrix(0,ncol=1,nrow=n)
    deta.dbeta <- matrix(0,ncol=k,nrow=n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dbeta0[i,]<-  B0[(i-m),] +  theta%*%deta.dbeta0[i-ma,]
      
      deta.dbeta[i,]<-  B[(i-m),] +  theta%*%deta.dbeta[i-ma,]
      
      deta.dphi[i,]<- A[(i-m),] +  theta%*%deta.dphi[i-ma,]
      
      deta.dtheta[i,]<- R[(i-m),] +  theta%*%deta.dtheta[i-ma,]
    }
    
    mM0 <- deta.dbeta0[(m+1):n,]
    mM <- deta.dbeta[(m+1):n,]
    pp <- deta.dphi[(m+1):n,]
    qq <- deta.dtheta[(m+1):n,]
    
    ymstar <- matrix((mustar),ncol=1)
    
    Ubeta0 <-  t(mM0) %*% mT %*% ymstar
    Ubeta <-  t(mM) %*% mT %*% ymstar
    Uphi <-    t(pp) %*% mT %*% ymstar
    Utheta <-  t(qq) %*% mT %*% ymstar
    
    ##############################################################################
    ##### START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    
    num1<-(mu-1)*exp(alpha)*log(1-tau)*critical.ly
    den1<- mu*(exp(alpha/mu)-exp(alpha))*(den.cr)^2
    num2.part1<-(y1-1)*critical.y
    num2.part2<-(log(1-tau)/den.cr-1)
    den2<-y1-y1*critical.y
    den3<- mu*(exp(alpha/mu)-exp(alpha))*den.cr
    Ualpha<-sum(1+1/alpha+num1/den1-num2.part1*num2.part2/den2+((mu-1)*exp(alpha))/den3-1/y1)
    
    ##### END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO alpha
    ##############################################################################
    
    return(c(Ubeta0,Ubeta,Uphi,Utheta,Ualpha))
  }#end score
  
  # initial values for estimation
  P.phi <- matrix(rep(NA,(n-m)*p1),ncol=p1) #ar
  for(i in 1:(n-m))
  {
    P.phi[i,] <- ynew[i+m-ar]
  }
  
  Q.theta <- matrix(rep(0,(n-m)*q1),ncol=q1) #ma
  
  x <- cbind(X[(m+1):n,],P.phi, Q.theta)#intercepto, covariáveis,ar, ma
  y1<-y[(m+1):n] 
  Ynew = linkfun(y1) 
  ajuste = lm.fit(x, Ynew) 
  
  mqo = c(ajuste$coef)
  for(i in 1:length(mqo))#porque no ma e MA vai virar NA, para então deixar 0
  {
    if (is.na(mqo[i])){
      mqo[i] <- 0
    }
  }
  library(GenSA)
  on.dmk.alpha<-function(alpha){-sum(log(dmk_alpha(alpha)))}
  
  gen.semchute <- GenSA(lower = c(.Machine$double.eps), 
                        upper = c(100),
                        fn = on.dmk.alpha, control=list(max.time=2))
  alpha<-gen.semchute$par
  reg <- c(mqo, alpha) # initializing the parameter values
  
  z <- c()
  opt.error<- tryCatch(optim(reg, loglik, score,
                             method = "BFGS",
                             control = list(fnscale = -1)), error = function(e) return("error"))  
  if(opt.error[1] == "error")
  {z$RMC=1
  #stop("optim error")
  warning("optim error")
  return(z)
  }
  opt <- optim(reg, loglik, score, 
               method = "BFGS", hessian = T, 
               control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
  #  opt2 <- optim(reg, loglik, method = "BFGS", hessian = T, control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
  #opt<-opt2
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE!")
    z$RMC=1
    return(z)
  }
  
  #  library(rootSolve)
  #print("verificando derivadas")
  #print(gradient(score,opt2$par))
  #print(hessian(loglik,opt2$par))
  
  # print(rbind(score(opt$par),gradient(loglik,opt$par)))
  # print(rbind(score(opt2$par),gradient(loglik,opt2$par)))
  # 
  z$conv <- opt$conv
  coef <- (opt$par)[1:(p1+q1+k+2)]
  z$coeff <- coef
  z$loglik <- opt$value
  beta0 <-coef[1]#intercept
  beta <-coef[2:(k+1)] #covariates
  phi<-coef[(k+2):(p1+k+1)] #ar
  theta <-coef[(p1+k+2):(p1+q1+k+1)] #ma
  alpha <-coef[(p1+q1+k+2)] # alpha parameter
  
  z$beta0 <- beta0
  z$beta <- beta
  z$phi <- phi
  z$theta <- theta
  z$alpha <- alpha
  
  z$ar=ar
  z$ma=ma
  z$AR=AR
  z$MA=MA
  z$delta=m
  z$roots=c(abs(polyroot(vector.root(z$ar,z$phi))),abs(polyroot(vector.root(z$ma,z$theta))))
  z$RMC=0
  if(any(z$roots<1)){warning("root(s) within the unity circle");z$RMC=1}
  
  errorhat <- rep(0,n) # E(error)=0
  etahat <- rep(NA,n)
  muhat0<- rep(NA,n)
  
  critical<-c()
  den.cr<-c()
  
  for(i in (m+1):n)
  {
    etahat[i] <- X[i,1]*beta0 + X[i,2:ncol(X)]%*%as.matrix(beta) + sum(phi*(ynew[i-ar]-X[i-ar,2:ncol(X)]%*%as.matrix(beta)) ) - sum(theta*errorhat[i-ma])
    errorhat[i] <- ynew[i]-etahat[i] #residuals on the predictor scale
  }
  
  muhat <- linkinv(etahat[(m+1):n])
  y1 <- y[(m+1):n]
  
  z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
  z$etahat <- etahat
  z$errorhat <- errorhat
  
  obs.inf<-function(y1,muhat)
  {
    critical<-c()
    den.cr<-c()
    for(i in 1:(n-m))
    {
      if(is.na(muhat[i])) muhat[i]<-.Machine$double.eps
      if(muhat[i]<.Machine$double.eps) muhat[i]<-.Machine$double.eps
      if(muhat[i]>0.9999999) muhat[i]<-0.9999999#1-.Machine$double.eps
      critical[i]<-(alpha-alpha/muhat[i])
      if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
      if(is.nan(critical[i])){critical[i]=-36.04365}
      if(critical[i] < -36.04365) critical[i] <- -36.04365
      den.cr[i]=log(1-exp(critical[i]))
      if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
    }
    
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU
    ###########################################################################################################
    e.ay=exp(alpha/y1)
    e.am=exp(alpha/muhat)
    critical.y=exp(alpha-alpha/y1)
    critical.ly<-log(1-critical.y)
    
    for (i in 1:length(y1)){
      if(is.infinite(e.ay[i])){e.ay[i]=.Machine$double.xmax}
      if(is.infinite(e.am[i])){e.am[i]=.Machine$double.xmax}
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
    }
    numerador<- ( exp(alpha)*alpha* ( den.cr*( (2*muhat*exp(alpha)+ e.am*(alpha-2*muhat))*log(1-tau)*critical.ly+exp(alpha)*alpha ) + (2*muhat*exp(alpha) +  e.am*(alpha-2*muhat))*den.cr*den.cr + 2*exp(alpha)*alpha*log(1-tau)*critical.ly ) ) 
    denominador <-( (muhat^4)*((exp(alpha)- e.am)^2) *den.cr*den.cr*den.cr )
    muhatstar.sec <-  numerador/denominador
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECTO TO MU   
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)# theta
    for(i in 1:(n-m))
    {
      for(j in 1:q1)
      {
        R[i,j] <- -sum(errorhat[i+m-ma[j]])
      }
    }
    
    A <- matrix(rep(NA,(n-m)*(p1)),ncol=(p1))# ar
    for(i in 1:(n-m))
    {
      for(j in 1:p1)
      {
        A[i,j] <- (ynew[i+m-ar[j]]-X[i+m-ar[j],2:ncol(X)]%*%as.matrix(beta)
        )
      }
    }    
    
    B0 <- matrix(rep(NA,(n-m)),ncol=1)#intercept
    for(i in 1:(n-m))
    {
      B0[i,1] <- X[i+m,1] 
    }
    
    B <- matrix(rep(NA,(n-m)*k),ncol=k)#covariates
    for(i in 1:(n-m))
    {
      for(j in 1:k)
      {
        B[i,j] <- X[i+m,1+j] - sum(phi*X[i+m-ar,1+j])
      }
    }
    
    bA <- array(NA,c(k,p1,(n-m)))# beta(ar)
    for(i in 1:(n-m))
    {
      for(b in 1:p1)
      {
        for(a in 1:k)
        {
          bA[a,b,i] <- -X[i+m-ar[b],(1+a)]
        }
      } 
    }
    
    deta.dbeta0 <- matrix(0,ncol=1,nrow=n)
    deta.dbeta <- matrix(0,ncol=k,nrow=n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta0beta0<- matrix(0, ncol=1,nrow=n)
    deta.dbeta0beta<- matrix(0, ncol=k,nrow=n)
    deta.dbeta0phi<- matrix(0, ncol=p1,nrow=n)
    deta.dbeta0theta<- matrix(0, ncol=q1,nrow=n)
    deta.dbetabeta<- array(0,dim=c(k,k,n))
    deta.dbetaphi <- array(0,dim=c(k,p1,n))
    deta.dbetatheta <- array(0,dim=c(k,q1,n))
    deta.dphiphi<-array(0,dim=c(p1,p1,n))
    deta.dphitheta<-array(0,dim=c(p1,q1,n))
    deta.dthetatheta<-array(0,dim=c(q1,q1,n))
    
    for(i in (m+1):n)
    {
      deta.dbeta0[i,]<- B0[(i-m),] +  theta%*%deta.dbeta0[i-ma,]
      deta.dbeta[i,]<- B[(i-m),] +  theta%*%deta.dbeta[i-ma,]
      deta.dphi[i,]<- A[(i-m),] +  theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] +  theta%*%deta.dtheta[i-ma,]
      deta.dbeta0beta0[i,]<- 0 +  theta%*%deta.dbeta0beta0[i-ma,]
      deta.dbeta0beta[i,]<-rep(0,k) +  theta%*%deta.dbeta0beta[i-ma,]
      deta.dbeta0phi[i,]<- rep(0,p1) +  theta%*%deta.dbeta0phi[i-ma,]
      for(a in 1:q1)
      {
        deta.dbeta0theta[i,a] = deta.dbeta0[i-ma[a],] +theta%*%deta.dbeta0theta[i-ma,a]
      }
      for(b in 1:k)
      {
        for(a in 1:k)
        {
          deta.dbetabeta[a,b,i] <- 0 + theta%*%deta.dbetabeta[a,b,i-ma]
        }
      }
      for(b in 1:p1)
      {
        for(a in 1:k)
        {
          deta.dbetaphi[a,b,i] <- bA[a,b,(i-m)] + theta%*%deta.dbetaphi[a,b,i-ma]
        }
      }
      for(b in 1:q1)
      {
        for(a in 1:k)
        {
          deta.dbetatheta[a,b,i] <- deta.dbeta[i-ma[b],a] + theta%*%deta.dbetatheta[a,b,i-ma]
        }
      }
      for(b in 1:p1)
      {
        for(a in 1:p1)
        {
          deta.dphiphi[a,b,i] <- 0 +  theta%*%deta.dphiphi[a,b,i-ma]
        }
      }
      for(b in 1:q1)
      {
        for(a in 1:p1)
        {
          deta.dphitheta[a,b,i]= deta.dphi[i-ma[b],a] +theta%*%deta.dphitheta[a,b,i-ma]
        }
      }
      for(b in 1:q1)
      {
        for(a in 1:q1)
        {
          deta.dthetatheta[a,b,i]= deta.dtheta[i-ma[a],b] +deta.dtheta[i-ma[b],a]+theta%*%deta.dthetatheta[a,b,i-ma]
        }
      }
    }
    
    mM0 <- matrix(deta.dbeta0[(m+1):n,],ncol=1,nrow=(n-m))
    pp <- matrix(deta.dphi[(m+1):n,], ncol=p1,nrow=(n-m))
    qq <- matrix(deta.dtheta[(m+1):n,], ncol=q1,nrow=(n-m))
    mM02<- deta.dbeta0beta0[(m+1):n,]
    mM02<- matrix(deta.dbeta0beta0[(m+1):n,], ncol=1,nrow=(n-m))
    B0p<- matrix(deta.dbeta0phi[(m+1):n,], ncol=p1,nrow=(n-m))
    B0q=matrix(deta.dbeta0theta[(m+1):n,], ncol=q1,nrow=(n-m))
    pp2<- array(deta.dphiphi[,,(m+1):n],dim=c(p1,p1,(n-m)))
    pq=array(deta.dphitheta[,,(m+1):n],dim=c(p1,q1,(n-m)))
    qq2=array(deta.dthetatheta[,,(m+1):n],dim=c(q1,q1,(n-m))) 
    mM <- matrix(deta.dbeta[(m+1):n,], ncol=k,nrow=(n-m))
    B0B<-matrix(deta.dbeta0beta[(m+1):n,], ncol=k,nrow=(n-m))
    mM2<- array(deta.dbetabeta[,,(m+1):n],dim=c(k,k,(n-m)))
    Bp<- array(deta.dbetaphi[,,(m+1):n],dim=c(k,p1,(n-m)))
    Bq=array(deta.dbetatheta[,,(m+1):n],dim=c(k,q1,(n-m)))
    ####START SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA (CONFERIDA, IGUAL AO ÚLTIMO TERMO DA HESSIANA se usar mu ao invés de muhat)
    ###########################################################################################################
    num1<-2*((muhat-1)^2)*exp(2*alpha)*log(1-tau)*critical.ly
    den1<-(muhat^2)*((exp(alpha)-e.am)^2)*((den.cr)^3)
    num2<-((muhat-1)^2)*exp(alpha)*log(1-tau)*critical.ly
    den2<-(muhat^2)*(e.am-exp(alpha))*((den.cr)^2)
    num3<-((muhat-1)^2)*exp(2*alpha)*log(1-tau)*critical.ly
    den3<-(muhat^2)*((exp(alpha)-e.am)^2)*((den.cr)^2)
    num4<-((muhat-1)^2)*exp(2*alpha)
    den4<-den3
    num5<-((muhat-1)^2)*exp(alpha)
    den5<-(muhat^2)*(e.am-exp(alpha))*den.cr
    num6<-num4
    den6<-(muhat^2)*((exp(alpha)-e.am)^2)*den.cr
    num7<-((y1-1)^2)*critical.y*(log(1-tau)/den.cr-1)
    den7<-(y1^2)*(1-critical.y)
    num8<-((y1-1)^2)*exp(2*alpha-2*alpha/y1)*(log(1-tau)/den.cr-1)
    den8<-(y1^2)*((critical.y-1)^2)
    num9<-2*(muhat-1)*exp(2*alpha)*(y1-1)*log(1-tau)
    den9<-muhat*y1*(exp(alpha)-e.am)*(e.ay-exp(alpha))*((den.cr)^2)
    
    Ualphaalpha<-sum(num1/den1+num2/den2+num3/den3+num4/den4+num5/den5+num6/den6-num7/den7-num8/den8+num9/den9-1/(alpha^2))
    ########################################################################################################### 
    ####END SECOND DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    ###########################################################################################################
    e.a=exp(alpha)
    numerador<-(e.a*(den.cr*(log(1-tau)*(muhat*(-e.a)*alpha*(y1-1)*(e.a-e.am)-y1*(muhat*e.a-e.am*(muhat*alpha+muhat-alpha))*(e.a-e.ay)*critical.ly)+(muhat-1)*e.a*alpha*y1*(e.a-e.ay))+2*(muhat-1)*e.a*alpha*y1*(e.a-e.ay)*log(1-tau)*critical.ly+y1*(muhat*e.a-e.am*(muhat*alpha+muhat-alpha))*(-(e.a-e.ay))*((den.cr) ^2)))
    denominador<-(muhat^3)*y1*((e.a-e.am)^2)*(e.a-e.ay)*((den.cr)^3)
    Umualpha<-numerador/denominador
    ########################################################################################################### 
    ####END DERIVATIVE FROM [DERIVATIVE FROM LOG LIKELIHOOD IN RESPECT TO MU] IN RESPECT TO ALPHA
    
    ####START DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU
    ###########################################################################################################
    critical.y=exp(alpha-alpha/y1)
    critical.ly<-log(1-critical.y)
    for (i in 1:length(y1)){
      if(is.infinite(critical.y[i])){critical.y[i]=.Machine$double.xmax}
      if(is.nan(critical.ly[i])){critical.ly[i]=-36.04365}
      if(critical.ly[i] < -36.04365) critical.ly[i] <- -36.04365 #para exp dar .Machine$double.eps
      
    }
    num1<-exp(alpha)*alpha*(den.cr+log(1-tau)*critical.ly)
    den1<-(muhat^2)*(exp(alpha/muhat)-exp(alpha))*((den.cr)^2)
    mustar<-num1/den1
    ########################################################################################################### 
    ####END DERIVATIVE FROM LOG LIKELIHOOD WITH RESPECT TO MU 
    
    mT <- diag(mu.eta(etahat[(m+1):n]))    
    num.mT2=dif2link(muhat)
    mT2=diag(-num.mT2*(mu.eta(etahat[(m+1):n])^3))
    mV <- diag(muhatstar.sec)
    mV0 <- diag(mustar)
    mualpha<-diag(Umualpha)
    vI <- matrix(rep(1,n-m),ncol=1)
    
    KB0B0 <- -(t(mM0)%*%mV%*%(mT^2)%*%mM0 + t(mM0)%*%mV0%*%mT2%*%mM0 + t(vI)%*%mV0%*%mT%*%mM02)
    KB0B <- -(t(mM0)%*%mV%*%(mT^2)%*%mM + t(mM0)%*%mV0%*%mT2%*%mM + t(vI)%*%mV0%*%mT%*%B0B)
    KB0p <- -(t(mM0)%*%mV%*%(mT^2)%*%pp + t(mM0)%*%mV0%*%mT2%*%pp + t(vI)%*%mV0%*%mT%*%B0p)
    KB0q <- -(t(mM0)%*%mV%*%(mT^2)%*%qq+t(mM0)%*%mV0%*%mT2%*%qq + t(vI)%*%mV0%*%mT%*%B0q)
    KB0alpha <- -t(mM0)%*% mualpha %*% mT %*% vI#
    KBB0<-t(KB0B)
    KBB=matrix(rep(NA,k*k),ncol=k)
    if(length(KBB)==1){
      KBB <- -(t(mM)%*%mV%*%(mT^2)%*%mM + t(mM)%*%mV0%*%mT2%*%mM + t(vI)%*%mV0%*%mT%*%mM2)
    }else{
      for(j in 1:k){
        for(i in 1:k){
          KBB[i,j] <- -(t(mM[,i])%*%mV%*%(mT^2)%*%mM[,j] + t(mM[,i])%*%mV0%*%mT2%*%mM[,j] + t(vI)%*%mV0%*%mT%*%as.matrix(mM2[i,j,]))
        }
      }
    }
    KBp=matrix(rep(NA,k*p1),ncol=p1)
    if(length(KBp)==1){
      KBp <- -(t(mM)%*%mV%*%(mT^2)%*%pp + t(mM)%*%mV0%*%mT2%*%pp + t(vI)%*%mV0%*%mT%*%Bp)
    }else{
      for(j in 1:p1){
        for(i in 1:k){
          
          KBp[i,j] <- -(t(mM[,i])%*%mV%*%(mT^2)%*%pp[,j] + t(mM[,i])%*%mV0%*%mT2%*%pp[,j] + t(vI)%*%mV0%*%mT%*%as.matrix(Bp[i,j,]))
        }
        
      }
    }
    KBq=matrix(rep(NA,k*q1),ncol=q1)
    if(length(KBq)==1){
      KBq <- -(t(mM)%*%mV%*%(mT^2)%*%qq+t(mM)%*%mV0%*%mT2%*%qq + t(vI)%*%mV0%*%mT%*%Bq)
    }else{
      for(j in 1:q1){
        for(i in 1:k){
          
          KBq[i,j] <- -(t(mM[,i])%*%mV%*%(mT^2)%*%qq[,j] + t(mM[,i])%*%mV0%*%mT2%*%qq[,j] + t(vI)%*%mV0%*%mT%*%as.matrix(Bq[i,j,]))
          
        }
      } 
    }
    KBalpha <- -t(mM)%*% mualpha %*% mT %*% vI
    
    KpB0 <- t(KB0p)
    KpB <- t(KBp)
    Kpp=matrix(rep(NA,p1*p1),ncol=p1)
    if(length(Kpp)==1){
      Kpp <- -(t(pp)%*%mV%*%(mT^2)%*%pp+t(pp)%*%mV0%*%mT2%*%pp + t(vI)%*%mV0%*%mT%*%pp2)
    }else{
      for(j in 1:p1){
        for(i in 1:p1){
          Kpp[i,j] <- -(t(as.matrix(pp[,i]))%*%mV%*%(mT^2)%*%as.matrix(pp[,j])+t(as.matrix(pp[,i]))%*%mV0%*%mT2%*%as.matrix(pp[,j]) + t(vI)%*%mV0%*%mT%*%as.matrix(pp2[i,j,]))
        }}
    }
    #print("Kpp");print(Kpp)
    Kpq=matrix(rep(NA,p1*q1),ncol=q1)
    #print("pq");print(pq)
    if(length(Kpq)==1){
      Kpq <- -(t(pp)%*%mV%*%(mT^2)%*%qq+t(pp)%*%mV0%*%mT2%*%qq + t(vI)%*%mV0%*%mT%*%pq)
    }else{
      for(j in 1:q1){
        for(i in 1:p1){
          
          Kpq[i,j] <- -(t(as.matrix(pp[,i]))%*%mV%*%(mT^2)%*%as.matrix(qq[,j])+t(as.matrix(pp[,i]))%*%mV0%*%mT2%*%as.matrix(qq[,j]) + t(vI)%*%mV0%*%mT%*%as.matrix(pq[i,j,]))
        }
      }
    }
    Kpalpha <- -t(pp) %*% mualpha %*% mT %*% vI
    
    KqB0 <- t(KB0q)
    KqB <- t(KBq)
    Kqp <- t(Kpq)
    Kqq=matrix(rep(NA,q1*q1),ncol=q1)
    if(length(Kqq)==1){
      Kqq <- -(t(qq)%*%mV%*%(mT^2)%*%qq+t(qq)%*%mV0%*%mT2%*%qq + t(vI)%*%mV0%*%mT%*%qq2)
    }else{
      for(j in 1:q1){
        for(i in 1:q1){
          Kqq[i,j] <- -(t(as.matrix(qq[,i]))%*%mV%*%(mT^2)%*%as.matrix(qq[,j])+t(as.matrix(qq[,i]))%*%mV0%*%mT2%*%as.matrix(qq[,j]) + t(vI)%*%mV0%*%mT%*%as.matrix(qq2[i,j,]))
        }}
    }
    Kqalpha <- -t(qq)%*% mualpha %*% mT %*% vI 
    
    KalphaB0 <- t(KB0alpha)
    KalphaB <- t(KBalpha)
    Kalphap <- t(Kpalpha)
    Kalphaq <- t(Kqalpha)
    Kalphaalpha <- -Ualphaalpha
    
    K <- rbind(
      cbind(KB0B0,KB0B,KB0p,KB0q,KB0alpha),
      cbind(KBB0,KBB,KBp,KBq,KBalpha),
      cbind(KpB0,KpB,Kpp,Kpq,Kpalpha),
      cbind(KqB0,KqB,Kqp,Kqq,Kqalpha),
      cbind(KalphaB0,KalphaB,Kalphap,Kalphaq,Kalphaalpha)
    )
    return(K)
  }
  K<-obs.inf(y1,muhat)
  # print("-hessiana = Matriz de informação observada condicional")
  # print(round(K,4))
  # print("hessiana numerica")
  # print(round(-opt$hessian,4))
  # print("comparando meu cálculo com hessiana da estimação numérica")
  # print(round((K+opt2$hessian),2))
  # print("comparando meu cálculo com hessiana numérica da estimação analítica")
  # print(round((K+opt$hessian),2))
  # print("soma diferença hessiana otimização numérica")
  # print(round(sum(abs(K+opt2$hessian)),2))
  # print("soma diferença hessiana numérica otimização analítica")
  # print(round(sum(abs(K+opt$hessian)),2))
  # 
  
  
  Ksolve<- tryCatch(solve(K), error = function(e) return("error"))
  
  if(Ksolve[1] == "error")
  {z$RMC=1#used at Monte-Carlo simulation for discard from the sample
  warning("Analytic Observed Information Matrix is not positive semi-definite")
  return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
  }else{sol=try(solve(K))}
  
  v<-diag(sol)#Variância assintótica dos esimadores
  
  for (i in 1:length(v))
  {
    if(is.na(v[i]) | is.nan(v[i]) | v[i]<0 )  {
      z$RMC=1
      warning("Analytic Observed Information Matrix is not positive semi-definite")
      return(z)#if Analytic Observed Information Matrix is not positive semi-definite, do not calculate
    }
  }
  #return(z) # # # # # # # # # # # # # # tirar depois, abaixo nao precisa na simulacao # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  
  z$zstat<-z$coeff/sqrt(v)
  resp<-rep(0,length(z$zstat))
  for (i in 1:length(resp)){
    if(abs(z$zstat[i])>qnorm(0.975))
    {
      resp[i] <- "H0 rejected"
    } else {resp[i] <- "H0 not rejected"}
  }
  #print(resp)
  #print("Intervalos de confiança")
  LI<-z$coeff-qnorm(0.975)*sqrt(v)
  LS<-z$coeff+qnorm(0.975)*sqrt(v)
  z$LI=LI
  z$LS=LS
  z$pvalues<-(1-pnorm(abs(z$zstat)))*2
  #print(pvalue)
  
  first_col<-c("intercept", 1:k,ar, ma,"alpha estimator")
  result <- matrix(c(first_col,round(c(z$coeff,z$zstat,LI,LS,z$pvalues),4),resp), nrow=length(z$coeff), ncol=7, byrow=F)
  colnames(result) <- c("Estimator","MLE","Wald's Statistic","Lower bound","Upper bound","p-value","Wald'S Test result")
  rownames(result)<-c("", rep("cov",k), rep("ar",length(ar)), rep("ma",length(ma)),"")
  
  #print(result,quot=F)
  z$coef.result<-result
  z$maic <- -2*(z$loglik)*(n/(n-m))+2*(length(opt$par)) 
  z$mbic <- -2*(z$loglik)*(n/(n-m))+(length(opt$par))*log(n)
  result2<-matrix(round(c(z$loglik,z$maic,z$mbic),4),nrow=3,ncol=1)
  rownames(result2)<-c("Log-likelihood","AIC","BIC")
  
  ###START error metrics
  mae<-sum(abs(y[(m+1):n]-z$fitted[(m+1):n]))/(n-m)
  
  sq<-rep(NA,n)
  for(i in (m+1):n)
  {
    sq[i]<-(y[i]-z$fitted[i])^2
  }  
  
  mse<-sum(sq[(m+1):n])/(n-m)
  
  rmse<-sqrt(mse)
  
  mape<-sum(abs((y[(m+1):n]-z$fitted[(m+1):n])/y[(m+1):n])*100)/(n-m)
  
  MdRAE<-median(abs(y[(m+1):n]-z$fitted[(m+1):n])/abs(y[(m+1):n]-y[(m+1-1):(n-1)]))#seasonal, if not, -1 not -S, If our model’s forecast equals to the benchmark’s forecast then the result is 1. If the benchmarks forecast are better than ours then the result will be above > 1. If ours is better than it’s below 1.
  MAEnaive<-sum(abs(y[(m+1+1):(n)]-y[(m+1):(n-1)]))/(n-m-1)#seasonal, if not, -1 not -S
  MASE<-mae/MAEnaive #Its value greater than one (1) indicates the algorithm is performing poorly compared to the naïve forecast.
  
  MAEnaive.star<-sum(abs(y[(m+1):(n)]-y[(m+1-1):(n-1)]))/(n-m)
  MASE.star<-mae/MAEnaive.star
  
  #Mean directional accuracy
  sign.y<-sign(y[(m+1):(n)]-y[(m):(n-1)])
  sign.f<-sign(z$fitted[(m+1):(n)]-y[(m):(n-1)])
  MDA.cont<-0
  for (i in (1):(n-m)){  
    if(sign.y[i]==sign.f[i]){MDA.cont<-MDA.cont+1}  
  }
  MDA<-MDA.cont/(n-m)
  
  MASE=MASE.star
  
  
  z$accuracyfitted<-accuracy<-matrix(round(c(mae,mse,rmse,mape,MdRAE,MASE,MDA),4), nrow=1, ncol=7, byrow=T)
  colnames(z$accuracyfitted) <-colnames(accuracy) <- c("MAE","MSE","RMSE","MAPE","MdRAE", "MASE","MDA")
  rownames(accuracy) <- c("")
  
  ytofit<-ts(c(y[1:n]),start=start(y),frequency=frequency(y))
  
  #in-sample 1:m fit
  #como fit
  etam <- c()
  errorm=errorhat
  for(i in m:1)
  {
    etam[i] <- X[i,1]*z$beta0 +X[i,2:ncol(X)]%*%as.matrix(z$beta)+ sum(phi*(ynew[i+ar])-X[i+ar,2:ncol(X)]%*%as.matrix(z$beta)) - sum(theta*errorm[i+ma])
    errorm[i] <- ynew[i]-etam[i] 
  }
  ym2 <- linkinv(etam[1:m])
  z$fitm2=ym2
  z$fit.all <- ts(c(ym2,muhat),start=start(y),frequency=frequency(y))
  
  if(steps!=0){
    #### out of sample forecast
    ynew_prev <- c(ynew,rep(NA,steps))
    y_prev <- c(z$fitted,rep(NA,steps))
    X_prev<-matrix(rep(1,(n+steps)), nrow=(n+steps), ncol=1, byrow=F)
    
    ntotal<-n+steps
    X_prev <- matrix(c(rep(1,ntotal),exvar), nrow=ntotal, ncol=(k+1),byrow = F)
    for(i in 1:steps)
    {
      ynew_prev[n+i] <- X_prev[n+i,1]*z$beta0 + X_prev[n+i,2:ncol(X_prev)]%*%as.matrix(z$beta) + sum(z$phi*(ynew_prev[n+i-ar]-X_prev[n+i-ar,2:ncol(X_prev)]%*%as.matrix(z$beta)
      ) ) - sum(z$theta*errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # residuals on the original scale y-mu 
    }
    
    z$forecast <- ts(c(rep(NA,n),y_prev[(n+1):(n+steps)]),start=start(y),frequency=frequency(y))
    
    fittedplusout_forecast <-  ts(c(rep(NA,m),muhat,y_prev[(n+1):(n+steps)]),start=start(y),frequency=frequency(y))
    #### rolling window forecast
    gy_prev <- c(rep(NA,n+steps))
    gy<-linkfun(y)
    yr_prev <- c(z$fitted,rep(NA,steps))
    
    for(i in 1:steps)
    {
      gy_prev[n+i] <- X_prev[n+i,1]*z$beta0 + X_prev[n+i,2:ncol(X_prev)]%*%as.matrix(z$beta) + sum(z$phi*(gy[n+i-ar]-X_prev[n+i-ar,2:ncol(X_prev)]%*%as.matrix(z$beta)
      ) ) - sum(z$theta*errorhat[n+i-ma])
      yr_prev[n+i] <- linkinv(gy_prev[n+i])
      errorhat[n+i] <- 0 # residuals on the original scale y-mu 
    }
    
    z$rollingforecast <- ts(c(rep(NA,n),yr_prev[(n+1):(n+steps)]),start=start(y),frequency=frequency(y))
    
  }
  
  z$serie <- y
  
  #quantile residuals
  
  critical<-c()
  den.cr<-c()
  for(i in (m+1):n)
  {
    if(is.na(z$fitted[i])) z$fitted[i]<-.Machine$double.eps
    if(z$fitted[i]<.Machine$double.eps) z$fitted[i]<-.Machine$double.eps
    if(z$fitted[i]>0.9999999) z$fitted[i]<-0.9999999#1-.Machine$double.eps
    critical[i]<-(alpha-alpha/z$fitted[i])
    if(is.na(critical[i])) critical[i]<- -.Machine$double.eps
    if(is.nan(critical[i])){critical[i]=-36.04365}
    if(critical[i] < -36.04365) critical[i] <- -36.04365
    den.cr[i]=log(1-exp(critical[i]))
    if(is.nan(den.cr[i])){den.cr[i]=-36.04365}
  }
  z$resid1 <- as.vector(qnorm(pmk(y[(m+1):n],alpha, log(1-tau)/den.cr[(m+1):n],log.p = FALSE ) ))
  #deviance residuals
  l_tilde <- (dmk(y[(m+1):n], alpha, log(1-tau)/log(1-exp(alpha-alpha/y[(m+1):n])), log = TRUE))#y[(m+1):n] where was mu
  l_hat <- (dmk(y[(m+1):n], alpha, log(1-tau)/den.cr[(m+1):n], log = TRUE))#z$fitted[(m+1):n] where was mu
  
  dt <- (l_tilde-l_hat)
  dt[which(dt<0)]<-0
  
  r2a<-sign(y[(m+1):n]-z$fitted[(m+1):n])
  r2b<-sqrt(2*(dt))
  z$resid2<-r2a*r2b
  
  z$deviance <- 2*sum(dt)
  z$dof.dev=(n-m-p1-q1)#desconsidera intercepto do eta e alpha da distribuição
  z$p_deviance <- 1 - pchisq(z$deviance, z$dof.dev)
  mresult<-matrix(round(c(z$loglik,z$maic,z$mbic,z$deviance),4),nrow=4,ncol=1)
  z$deviance.star <- 2*sum(dt)*n/(n-m)
  z$dof.dev=(n-m-length(opt$par)-1)
  rownames(mresult)<-c("Log-likelihood","AIC","BIC","Deviance")
  colnames(mresult)<-c("")
  z$mresult<-mresult
  
  if(resid==1) {
    residual <- z$resid1
    
  }
  
  if(resid==2) {
    residual <- z$resid2
  }
  
  z$residual<-residual
  ###################################################
  ######### GRAPHICS ################################
  
  if(graph==T)
  {
    t<-seq(-5,n+6,by=1)
    w1<-5
    h1<-4
    pdf("resid_v_ind.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1))
      par(mgp=c(1.7, 0.45, 0))
      plot(residual,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
      lines(t,rep(-3,n+12)#length(residual))
            ,lty=2,col=1)
      lines(t,rep(3,n+12)#length(residual))
            ,lty=2,col=1)
      lines(t,rep(-2,n+12)#length(residual))
            ,lty=3,col=1)
      lines(t,rep(2,n+12)#length(residual))
            ,lty=3,col=1)
    }
    dev.off()
    
    pdf("resid_v_fitted.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      plot(as.vector(z$fitted[(m+1):n]),as.vector(residual), main=" ", pch = "+",
           xlab="Fitted values",ylab="Residuals",ylim=c(-4,4))
      lines(t,rep(-3,n+12),lty=2,col=1)
      lines(t,rep(3,n+12),lty=2,col=1)
      lines(t,rep(-2,n+12),lty=3,col=1)
      lines(t,rep(2,n+12),lty=3,col=1)
    }
    dev.off()
    pdf("obs_v_fit.pdf",width=5, height=4)### abre no navegador google chrome só
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      plot(as.vector(z$fitted), as.vector(ytofit), main=" ", pch = "+",
           xlab="Fitted values",ylab="Observed data",
           xlim=c(0.95*min(y),max(y)*1.05),
           ylim=c(0.95*min(y),max(y)*1.05))
      lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    }
    dev.off()
    
    pdf("resid_density.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(1.5, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      densidade<-density(residual)
      plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
      lines(densidade$x,dnorm(densidade$x),lty=2)
      legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
             pt.bg="white", lty=c(1,2), bty="n")
    }
    dev.off()
    pdf("resid_FAC.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      acf(residual,ylab="ACF",xlab="Lag") 
    }
    dev.off()
    pdf("resid_FACP.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      pacf(residual,ylab="PACF",xlab="Lag")
    }
    dev.off()
    pdf("qq_plot.pdf",width=5, height=4)
    {  
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) 
      par(mgp=c(1.7, 0.45, 0))
      qqnorm(residual, pch = "+",
             xlim=c(0.95*min(residual),max(residual)*1.05),
             ylim=c(0.95*min(residual),max(residual)*1.05),
             main="",xlab="Normal quantiles",ylab="Empirical quantiles")
      lines(c(-10,10),c(-10,10),lty=2)
    }
    dev.off()
    pdf("adjusted.pdf",width=5, height=4)
    {
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
      par(mgp=c(1.7, 0.45, 0))
      plot(ytofit,type="l",ylab="Serie",xlab="Time")
      lines(z$fitted,col="blue",lty=2)
      legend("bottomleft",c("Observed data","Fitted values"),#pch=vpch,
             pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
    }
    dev.off()
    pdf("fitted.all.pdf",width=5, height=4)
    {
      sfit=start(y)[1]+start(y)[2]/12+(m-1)/12
      par(mfrow=c(1,1))
      par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
      par(mgp=c(1.7, 0.45, 0))
      plot(ytofit,type="l",ylab="Serie",xlab="Time")
      abline(v=sfit,lty=2)
      lines(z$fit.all,col="blue",lty=2)
      legend("bottomleft",c("Observed data","Fitted values"),#pch=vpch,
             pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
    }
    dev.off()
    if(steps!=0){
      pdf("fittedforecast.pdf",width=5, height=4)
      {
        fim<-end(y)[1]+end(y)[2]/12
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev,type="l",col="blue",lty=2, ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
        abline(v=fim,lty=2)
        abline(v=n,lty=2)
        lines(as.vector(y))
        legend("bottomleft",c("Observed data","Fitted and forecast values"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
      }
      dev.off()
      
      pdf("forecast.pdf",width=5, height=4)
      {
        fim<-end(y)[1]+end(y)[2]/12
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev[(n+1):(n+steps)],type="l",col="blue",lty=2, ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
        abline(v=fim,lty=2)
        abline(v=n,lty=2)
        lines(as.vector(y[(n+1):(n+steps)]))
        legend("bottomleft",c("Observed data","Forecast values"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
      }
      dev.off()
    }
  }#END GRAPHICS
  
  
  ########################################################################
  ########################   residual analysis   ########################
  ########################################################################
  
  #null hypothesis: non-autocorrelation
  
  if(AR[1]!=0 | MA[1]!=0){
    aclag <- 2*S
  }
  
  ljungbox<- Box.test(residual, lag = aclag, type = "Ljung-Box", fitdf=(p1+q1))
  boxpierce<-Box.test(residual, lag = aclag, type = c( "Box-Pierce"), fitdf=(p1+q1))#Box-Pierce test
  monti<- Monti.test(residual, lag = aclag, type = "Ljung-Box", fitdf=(p1+q1))#Ljung-Box test as the Box.test function using pacf replacing acf
  
  z$boxpierce<-boxpierce$statistic
  z$p_boxpierce<-boxpierce$p.value
  
  z$ljungbox<-ljungbox$statistic
  z$p_ljungbox<-ljungbox$p.value
  
  z$monti<-monti$statistic
  z$p_monti<-monti$p.value
  
  #null hypothesis: normality
  library(tseries)
  jarquebera<-jarque.bera.test(residual)
  z$jarquebera<-jarquebera$statistic
  z$p_jarquebera<-jarquebera$p.value
  
  library(nortest)
  andersondarling=ad.test(residual)
  z$andersondarling<-andersondarling$statistic
  z$p_andersondarling<-andersondarling$p.value
  
  #null hypothesis: non-heteroscedasticity (constant variance)
  
  library(FinTS)
  arch<-ArchTest(residual, lags=10) 
  z$arch<-arch$statistic
  z$p_arch<-arch$p.value
  
  #null hypothesis: no stable seasonality
  library(seastests)
  friedman=fried(ts(y[1:n],frequency=S))
  z$friedman=friedman$stat
  z$p_friedman=friedman$Pval
  #if the null hipothesis is rejected at the 1% significance level then the series is considered to be seasonal
  
  ########################################################################
  ########################   forecast analysis   ########################
  ########################################################################
  if(steps!=0){
    if(validation==T){
      ###START FORECAST error metrics
      accuracyforecast<-function(y_prev,steps){
        maef<-sum(abs(y[(n+1):(n+steps)]-y_prev[(n+1):(n+steps)]))/(steps)
        
        sqf<-rep(NA,steps)
        
        for(i in 1:steps)
        {
          sqf[i]<-(y[n+i]-y_prev[n+i])^2
        }  
        
        #print((y[(n+1):(n+steps)]-y_prev[(n+1):(n+steps)])^2)
        msef<-sum(sqf)/steps
        
        rmsef<-sqrt(msef)
        
        mapef<-sum(abs((y[(n+1):(n+steps)]-y_prev[(n+1):(n+steps)])/y[(n+1):(n+steps)])*100)/steps
        
        MdRAEf<-median(abs(y[(n+1):(n+steps)]-y_prev[(n+1):(n+steps)])/abs(y[(n+1):(n+steps)]-y[(n+1-S):(n+steps-S)]))#seasonal, if not, -1 not -S, If our model’s forecast equals to the benchmark’s forecast then the result is 1. If the benchmarks forecast are better than ours then the result will be above > 1. If ours is better than it’s below 1.
        
        MAEnaivef<-sum(abs(y[(n+S+1):(n+steps)]-y[(n+1):(n+steps-S)]))/(steps-S)#seasonal, if not, -1 not -S
        
        MASEf<-maef/MAEnaivef #Its value greater than one (1) indicates the algorithm is performing poorly compared to the naïve forecast.
        
        MAEnaivef.star<-sum(abs(y[(n+1):(n+steps)]-y[(n+1-S):(n+steps-S)]))/steps
        MASEf.star<-maef/MAEnaivef.star
        
        #Mean directional accuracy
        sign.yf<-sign(y[(n+1):(n+steps)]-y[(n):(n+steps-1)])
        sign.ff<-sign(y_prev[(n+1):(n+steps)]-y[(n):(n+steps-1)])
        MDAf.cont<-0
        for (i in 1:steps){   
          if(sign.yf[i]==sign.ff[i]){MDAf.cont<-MDAf.cont+1}  
        }
        
        MDAf<-MDAf.cont/steps
        
        MASEf=MASEf.star
        accuracyf<-matrix(round(c(maef,msef,rmsef,mapef,MdRAEf,MASEf,MDAf),4), nrow=1, ncol=7, byrow=T)
        colnames(accuracyf) <- c("MAE","MSE","RMSE","MAPE","MdRAE","MASE","MDA")
        rownames(accuracyf) <- c("Accuracy forecast")
        return(accuracyf)
      }
      accuracytraditionalforecast<-accuracyrollingwindow<-matrix(rep(NA,7*steps),nrow=steps, ncol=7, byrow=T)
      colnames(accuracytraditionalforecast) <- colnames(accuracyrollingwindow) <- c("MAE","MSE","RMSE","MAPE","MdRAE","MASE","MDA")
      rownames(accuracytraditionalforecast) <- rownames(accuracyrollingwindow) <- 1:steps
      for (i in 1:steps){
        accuracytraditionalforecast[i,]<-accuracyforecast(y_prev,steps=i)
        accuracyrollingwindow[i,]<-accuracyforecast(yr_prev,steps=i)
      }
      z$accuracyforecast<-accuracytraditionalforecast
      z$accuracyrollingwindow<-accuracyrollingwindow
    }
  }
  
  rownames(z$accuracyfitted) <-rownames(accuracy) <- c("Fitted accuracy")
  
  diagnostic<-matrix(round(c(z$boxpierce,z$ljungbox,z$monti,z$jarquebera,z$andersondarling,z$arch,
                             z$p_boxpierce,z$p_ljungbox,z$p_monti,z$p_jarquebera,z$p_andersondarling,z$p_arch
  ),4), nrow=2, ncol=6, byrow=T)
  colnames(diagnostic) <- c("Box-Pierce test","Ljung-Box tes","Monti test","Jarque-Bera test","Anderson-Darling test","Arch test")
  rownames(diagnostic) <- c("Statistic","P-value")
  z$diagnostic=diagnostic
  if(print==T){
    print("MKSARMAX",quote=F)
    print(z$coef.result,quote=F)
    message("")
    print(c("Log-likelihood =",round(z$loglik,4)),quote=F)
    print(c("MAIC =",round(z$maic,4),"MBIC =",round(z$mbic,4)),quote=F)
    print(c("Deviance =",round(z$deviance,4)," DF:",z$dof.dev,"Deviance* =",round(z$deviance.star,4)),quote=F)
    message("")  
    if(resid==1) {
      print("Quantile residuals:",quote=F)
    }
    
    if(resid==2) {
      print("Deviance residuals:",quote=F)
    }
    print(summary(z$residual))
    message("")
    print(z$diagnostic)
    message("")
    print(z$accuracyfitted)
    message("")
    if(steps!=0 & validation==T){
      print("Traditional forecast accuracy:",quote=F)
      print(z$accuracyforecast)
      message("")
      print("Rolling window forecast accuracy:",quote=F)
      print(z$accuracyrollingwindow)
    }
  }
  if(check==TRUE){
    opt2 <- optim(reg, loglik, method = "BFGS", hessian = T, control = list(fnscale = -1))    
    library(rootSolve)
    print("verificando derivadas")
    print(gradient(score,opt2$par))
    print(hessian(loglik,opt2$par))
    print("rbind(score(opt$par),gradient(loglik,opt$par))")
    print(rbind(score(opt$par),gradient(loglik,opt$par)))
    print("rbind(score(opt2$par),gradient(loglik,opt2$par))")
    print(rbind(score(opt2$par),gradient(loglik,opt2$par)))
    print("-hessiana = Matriz de informação observada condicional")
    print(round(K,4))
    print("hessiana numerica")
    print(round(-opt$hessian,4))
    print("comparando meu cálculo com hessiana da estimação numérica")
    print(round((K+opt2$hessian),2))
    print("comparando meu cálculo com hessiana numérica da estimação analítica")
    print(round((K+opt$hessian),2))
    print("soma diferença hessiana otimização numérica")
    print(round(sum(abs(K+opt2$hessian)),2))
    print("soma diferença hessiana numérica otimização analítica")
    print(round(sum(abs(K+opt$hessian)),2))
  }
  return(z)
}#fim estimação
