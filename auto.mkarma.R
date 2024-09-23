#max.order: maximum order p,q,P,Q to check first model
#till 4

auto.mkarma<-function(y,steps=10,exvar=matrix(NA, nrow=1, ncol=1, byrow=F),max.order=1,graph=F,print=F,type="Box-Jenkins",sort="mbic",S=12,resid=1){
  library(tictoc)
  
  nome1<-paste("auto.mkarma.txt", sep="")
  
  tic("time")
  sink(nome1)
  library(seastests)
  friedman=fried(ts(y[1:(length(y)-steps)],frequency=S))#Friedman Rank test
  source("fit.R")
  set.seed(1934)
  ar=c(1)
  ma=c(0)
  AR=c(0)
  MA=c(0)
  m.tenth=m.ninth=m.eighth=m.seventh=m.sixth=m.fifth=m.fourth=m.third=m.second=m.first<-c()
  first<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
  
  print(first$RMC)
  tenth=ninth=eighth=seventh=sixth=fifth=fourth=third=second=c()
  tenth$ar=ninth$ar=eighth$ar=seventh$ar=sixth$ar=fifth$ar=fourth$ar=third$ar=second$ar=first$ar<-ar
  tenth$ma=ninth$ma=eighth$ma=seventh$ma=sixth$ma=fifth$ma=fourth$ma=third$ma=second$ma=first$ma<-ma
  tenth$AR=ninth$AR=eighth$AR=seventh$AR=sixth$AR=fifth$AR=fourth$AR=third$AR=second$AR=first$AR<-AR
  tenth$MA=ninth$MA=eighth$MA=seventh$MA=sixth$MA=fifth$MA=fourth$MA=third$MA=second$MA=first$MA<-MA
  tenth$mbic=ninth$mbic=eighth$mbic=seventh$mbic=sixth$mbic=fifth$mbic=fourth$mbic=second$mbic=third$mbic<-100
  tenth$accuracyfitted=ninth$accuracyfitted=eighth$accuracyfitted=seventh$accuracyfitted=sixth$accuracyfitted=fifth$accuracyfitted=fourth$accuracyfitted=second$accuracyfitted=third$accuracyfitted<-matrix(rep(100,7),ncol=7)
  if(first$RMC==1){
    first$mbic<-100
    first$accuracyfitted<-matrix(rep(100,7),ncol=7)
  }
  if(type=="Box-Jenkins"){
    #start seasonal models
    if (friedman$Pval<=0.01 ){ 
      print(paste("Friedman Rank test. P-value=",friedman$Pval))
      print(paste("Including seasonal models on search for the best model."))
      n <- 4
      if (max.order==4){
        l <- rep(list(0:4), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]
         
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if(ar==3){ar=c(1,2,3)}else if(ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)}else if(ma==3){ma=c(1,2,3)}else if(ma==4){ma=c(1,2,3,4)}
          if (AR==2){AR=c(1,2)}else if(AR==3){AR=c(1,2,3)}else if(AR==4){AR=c(1,2,3,4)}
          if (MA==2){MA=c(1,2)}else if(MA==3){MA=c(1,2,3)}else if(MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 & all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }
        }
        
        print("first model up to order 4");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second first model up to order 4")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third first model up to order 4")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth first model up to order 4")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth first model up to order 4")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth first model up to order 4")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh first model up to order 4")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth first model up to order 4")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth first model up to order 4")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth first model up to order 4")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
      }
      
      if (max.order==3){
        l <- rep(list(0:3), n)
        possible<-expand.grid(l)
        
        models1<-possible[3:nrow(possible),]
         
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,2,3)}
          if (AR==2){AR=c(1,2)}else if (AR==3){AR=c(1,2,3)}
          if (MA==2){MA=c(1,2)}else if (MA==3){MA=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }
        }
        
        print("first model up to order 3");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second first model up to order 3")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third first model up to order 3")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth first model up to order 3")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth first model up to order 3")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth first model up to order 3")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh first model up to order 3")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth first model up to order 3")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth first model up to order 3")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth first model up to order 3")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
      }
      
      if (max.order==2){
        l <- rep(list(0:2), n)
        possible<-expand.grid(l)
        
        models1<-possible[3:nrow(possible),]
         
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          if (AR==2){AR=c(1,2)}
          if (MA==2){MA=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }
        }
        
        
        print("first model up to order 2");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second first model up to order 2")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third first model up to order 2")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth first model up to order 2")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth first model up to order 2")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth first model up to order 2")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh first model up to order 2")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth first model up to order 2")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth first model up to order 2")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth first model up to order 2")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
      }
      
      if (max.order==1){
        l <- rep(list(0:1), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]
         
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }
        }
        print("first model up to order 1");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second first model up to order 1")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third first model up to order 1")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth first model up to order 1")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth first model up to order 1")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth first model up to order 1")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh first model up to order 1")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth first model up to order 1")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth first model up to order 1")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth first model up to order 1")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
      }
      
      #end seasonal models
    }else{#start non-seasonal models
      print(paste("Friedman Rank test. P-value=",friedman$Pval))
      print(paste("Excluding seasonal models on search for the best model."))
      n <- 2
      if (max.order==4){
        l <- rep(list(0:4), n)
        possible<-expand.grid(l)
        
        models1<-possible[3:nrow(possible),]
 
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,2,3)}else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,2,3)}else if (ma==4){ma=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }
        }
        
        print("first model up to order 4");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        print("second model up to order 4")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        print("third model up to order 4")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        print("fourth model up to order 4")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        print("fifth model up to order 4")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        print("sixth model up to order 4")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        print("seventh model up to order 4")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        print("eighth model up to order 4")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        print("ninth model up to order 4")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        print("tenth model up to order 4")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
      }
      
      if (max.order==3){
        l <- rep(list(0:3), n)
        possible<-expand.grid(l)
        
        models1<-possible[3:nrow(possible),]
        
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }
        }
        
        print("first model up to order 3");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        print("second model up to order 3")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        print("third model up to order 3")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        print("fourth model up to order 3")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        print("fifth model up to order 3")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        print("sixth model up to order 3")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        print("seventh model up to order 3")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        print("eighth model up to order 3")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        print("ninth model up to order 3")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        print("tenth model up to order 3")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
      }
      
      if (max.order==2){
        l <- rep(list(0:2), n)
        possible<-expand.grid(l)
        
        models1<-possible[3:nrow(possible),]
        
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=")
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }
        }
        
        print("first model up to order 2");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        
        print("second model up to order 2")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        
        print("third model up to order 2")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        
        print("fourth model up to order 2")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        
        print("fifth model up to order 2")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        
        print("sixth model up to order 2")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        
        print("seventh model up to order 2")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        
        print("eighth model up to order 2")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        
        print("ninth model up to order 2")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        
        print("tenth model up to order 2")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        
      }
      
      if (max.order==1){
        l <- rep(list(0:1), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]
      
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }
        }
        print("first model up to order 1");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        
        print("second model up to order 1")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        
        print("third model up to order 1")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        
        print("fourth model up to order 1")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        
        print("fifth model up to order 1")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        
        print("sixth model up to order 1")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        
        print("seventh model up to order 1")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        
        print("eighth model up to order 1")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        
        print("nint model up to order 1")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        
        print("tenth model up to order 1")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        
      }
    }
  }else if(type=="all"){
    #start seasonal models
    if (friedman$Pval<=0.01 ){ 
      print(paste("Friedman Rank test. P-value=",friedman$Pval))
      
      print(paste("Including seasonal models on search for the best model."))
      n <- 4
      
      if (max.order==4){
        
        l <- rep(list(0:4), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]
        for (i in 1:nrow(models1)){#first step: 1-2-3-4
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second=first
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third=second
              second<-mod
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth=third
              third<-mod
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth=fourth
              fourth<-mod
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth=fifth
              fifth<-mod
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh=sixth
              sixth<-mod
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth
              ninth=eighth
              eighth=seventh
              seventh<-mod
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth
              ninth=eighth
              eighth<-mod
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth
              ninth<-mod
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        sep2<-which(possible[,1]==2 | possible[,2]==2 | possible[,3]==2 | possible[,4]==2)
        models2<-possible[sep2,]
          
        for (i in 1:nrow(models2)){#first step: 1-1,2-3-4
          order<-models2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          if (AR==2){AR=c(1,2)}
          if (MA==2){MA=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        sep3<-which(possible[,1]==3 | possible[,2]==3 | possible[,3]==3 | possible[,4]==3)
        models3<-possible[sep3,]
        
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,2,3-4
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,2,3)}
          if (AR==2){AR=c(1,2)}else if (AR==3){AR=c(1,2,3)}
          if (MA==2){MA=c(1,2)}else if (MA==3){MA=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,3-4
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,3)}
          if (AR==2){AR=c(1,2)}else if (AR==3){AR=c(1,3)}
          if (MA==2){MA=c(1,2)}else if (MA==3){MA=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        for (i in 1:nrow(models3)){#first step: 1-1,2-2,3-4
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(2,3)}
          if (AR==2){AR=c(1,2)}else if (AR==3){AR=c(2,3)}
          if (MA==2){MA=c(1,2)}else if (MA==3){MA=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        sep3.2=which(models3[,1]==2 | models3[,2]==2 | models3[,3]==2 | models3[,4]==2)
        if(length(sep3.2)!=0){
        models3.2=models3[sep3.2,]      
        
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,3-4
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)}
          if (ma==3){ma=c(1,3)}
          if (AR==3){AR=c(1,3)}
          if (MA==3){MA=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,2,3-4
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)}
          if (ma==3){ma=c(1,2,3)}
          if (AR==3){AR=c(1,2,3)}
          if (MA==3){MA=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        for (i in 1:nrow(models3.2)){#first step: 1-2-2,3-4
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)}
          if (ma==3){ma=c(2,3)}
          if (AR==3){AR=c(2,3)}
          if (MA==3){MA=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }}
          }}
        }
        
        l2 <- rep(list(0:2), n)
        possible2<-expand.grid(l2)
        for (i in 1:nrow(possible2)){
          order<-possible2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,2)}
            if(ma==1){ma=c(1,2)}
            if(AR==1){AR=c(1,2)}
            if(MA==1){MA=c(1,2)}
             mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        l3 <- rep(list(0:4), n)
        possible3<-expand.grid(l3)
        for (i in 1:nrow(possible3)){
          order<-possible3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,3)}else if(ar==2){ar=c(2,3)} else if(ar==4){ar=c(1,2,3)}
            if(ma==1){ma=c(1,3)}else if(ma==2){ma=c(2,3)} else if(ma==4){ma=c(1,2,3)}
            if(AR==1){AR=c(1,3)}else if(AR==2){AR=c(2,3)} else if(AR==4){AR=c(1,2,3)}
            if(MA==1){MA=c(1,3)}else if(MA==2){MA=c(2,3)} else if(MA==4){MA=c(1,2,3)}
             mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        l4 <- rep(list(0:8), n)
        possible4<-expand.grid(l4)
        for (i in 1:nrow(possible4)){
          order<-possible4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,4)}else if(ar==2){ar=c(2,4)}else if(ar==3){ar=c(3,4)}else if(ar==5){ar=c(1,2,4)}else if(ar==6){ar=c(1,3,4)}else if(ar==7){ar=c(2,3,4)}else if(ar==8){ar=c(1,2,3,4)}
            if(ma==1){ma=c(1,4)}else if(ma==2){ma=c(2,4)}else if(ma==3){ma=c(3,4)}else if(ma==5){ma=c(1,2,4)}else if(ma==6){ma=c(1,3,4)}else if(ma==7){ma=c(2,3,4)}else if(ma==8){ma=c(1,2,3,4)}
            if(AR==1){AR=c(1,4)}else if(AR==2){AR=c(2,4)}else if(AR==3){AR=c(3,4)}else if(AR==5){AR=c(1,2,4)}else if(AR==6){AR=c(1,3,4)}else if(AR==7){AR=c(2,3,4)}else if(AR==8){AR=c(1,2,3,4)}
            if(MA==1){MA=c(1,4)}else if(MA==2){MA=c(2,4)}else if(MA==3){MA=c(3,4)}else if(MA==5){MA=c(1,2,4)}else if(MA==6){MA=c(1,3,4)}else if(MA==7){MA=c(2,3,4)}else if(MA==8){MA=c(1,2,3,4)}
             mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA=");print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        sep4<-which(possible[,1]==4 | possible[,2]==4| possible[,3]==4| possible[,4]==4)
        models4<-possible[sep4,]
  
        for (i in 1:nrow(models4)){#first step 1-2-3-1,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==4){ar=c(1,4)}
          if (ma==4){ma=c(1,4)}
          if (AR==4){AR=c(1,4)}
          if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-3-2,4
        for (i in 1:nrow(models4)){#first step 1-2-3-2,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==4){ar=c(2,4)}
          if (ma==4){ma=c(2,4)}
          if (AR==4){AR=c(2,4)}
          if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-3-3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==4){ar=c(3,4)}
          if (ma==4){ma=c(3,4)}
          if (AR==4){AR=c(3,4)}
          if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-3-1,2,4
        for (i in 1:nrow(models4)){#first step 1-2-3-1,2,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==4){ar=c(1,2,4)}
          if (ma==4){ma=c(1,2,4)}
          if (AR==4){AR=c(1,2,4)}
          if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-3-1,3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-1,3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==4){ar=c(1,3,4)}
          if (ma==4){ma=c(1,3,4)}
          if (AR==4){AR=c(1,3,4)}
          if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-3-1,2,3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-1,2,3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==4){ar=c(1,2,3,4)}
          if (ma==4){ma=c(1,2,3,4)}
          if (AR==4){AR=c(1,2,3,4)}
          if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-3-2,3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-2,3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==4){ar=c(2,3,4)}
          if (ma==4){ma=c(2,3,4)}
          if (AR==4){AR=c(2,3,4)}
          if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        sep4.2=which(models4[,1]==2 | models4[,2]==2 | models4[,3]==2 | models4[,4]==2)
        if(length(sep4.2)!=0){
        models4.2=models4[sep4.2,]
        
        #first step 1-1,2-3-1,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,4)}
          if (AR==2){AR=c(1,2)} else if (AR==4){AR=c(1,4)}
          if (MA==2){MA=c(1,2)} else if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-3-2,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-2,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==4){AR=c(2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-3-3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==4){AR=c(3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-3-1,2,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,2,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==4){AR=c(1,2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-3-1,3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==4){AR=c(1,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-3-1,2,3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,2,3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==4){AR=c(1,2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-3-2,3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-2,3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==4){AR=c(2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        }
        sep4.3=which(models4[,1]==3 | models4[,2]==3 | models4[,3]==3 | models4[,4]==3)
        if(length(sep4.3)!=0){
        models4.3=models4[sep4.3,]
        
        #first step 1-2-1,2,3-1,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,4)}
          if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,4)}
          if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,2,3-2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,4)}
          if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(2,4)}
          if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,2,3-3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(3,4)}
          if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(3,4)}
          if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,2,3-1,2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,4)}
          if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,2,4)}
          if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,2,3-1,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,3,4)}
          if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,3,4)}
          if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,2,3-1,2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,3,4)}
          if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,2,3,4)}
          if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,2,3-2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,3,4)}
          if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(2,3,4)}
          if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        #first step 1-2-1,3-1,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,4)}
          if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,4)}
          if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,3-2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,4)}
          if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(2,4)}
          if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,3-3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(3,4)}
          if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(3,4)}
          if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,3-1,2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,4)}
          if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,2,4)}
          if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,3-1,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,3,4)}
          if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,3,4)}
          if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,3-1,2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,3,4)}
          if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,2,3,4)}
          if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-1,3-2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,3,4)}
          if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(2,3,4)}
          if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-2,3-1,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,4)}
          if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,4)}
          if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-2,3-2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,4)}
          if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(2,4)}
          if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-2,3-3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(3,4)}
          if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(3,4)}
          if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-2,3-1,2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,4)}
          if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,2,4)}
          if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-2,3-1,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,3,4)}
          if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,3,4)}
          if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-2,3-1,2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,3,4)}
          if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,2,3,4)}
          if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-2-2,3-2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,3,4)}
          if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(2,3,4)}
          if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        sep4.3.2=which(models4.3[,1]==2 | models4.3[,2]==2 | models4.3[,3]==2 | models4.3[,4]==2)
        if(length(sep4.3.2)!=0){
        models4.3.2=models4.3[sep4.3.2,]
        
        #first step 1-1,2-1,2,3-1,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-1,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,2,3-2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,2,3-3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,2,3-1,2,4
        for (i in 1:nrow(models4.3.2)){#first step  1-1,2-1,2,3-1,2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,2,3-1,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-1,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,2,3-1,2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-1,2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(1,2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,2,3-2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,2,3)} else if (AR==4){AR=c(2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,2,3)} else if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,3-1,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,3-2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,3-3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,3-1,2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,3-1,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,3-1,2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(1,2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-1,3-2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(1,3)} else if (AR==4){AR=c(2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(1,3)} else if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-2,3-1,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-2,3-2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-2,3-3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-2,3-1,2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,2,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-2,3-1,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-2,3-1,2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(1,2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        #first step 1-1,2-2,3-2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,3,4)}
          if (AR==2){AR=c(1,2)} else if (AR==3){AR=c(2,3)} else if (AR==4){AR=c(2,3,4)}
          if (MA==2){MA=c(1,2)} else if (MA==3){MA=c(2,3)} else if (MA==4){MA=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        }
        }
        print("first model up to order 4");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second model up to order 4")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third model up to order 4")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth model up to order 4")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth model up to order 4")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth model up to order 4")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh model up to order 4")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth model up to order 4")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth model up to order 4")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth model up to order 4")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
      }
      
      if (max.order==3){
        
        l <- rep(list(0:3), n)
        possible<-expand.grid(l)#255 nmodels
        models1<-possible[3:nrow(possible),]
         for (i in 1:nrow(models1)){#first step: 1-2-3
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        sep2<-which(possible[,1]==2 | possible[,2]==2 | possible[,3]==2 | possible[,4]==2)
        models2<-possible[sep2,]#175 models
       
        for (i in 1:nrow(models2)){#first step: 1-1,2-3
          order<-models2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          if (AR==2){AR=c(1,2)}
          if (MA==2){MA=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        sep3<-which(possible[,1]==3 | possible[,2]==3 | possible[,3]==3 | possible[,4]==3)
        models3<-possible[sep3,]
         
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,2,3
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,2,3)}
          if (AR==2){AR=c(1,2)}else if (AR==3){AR=c(1,2,3)}
          if (MA==2){MA=c(1,2)}else if (MA==3){MA=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,3
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,3)}
          if (AR==2){AR=c(1,2)}else if (AR==3){AR=c(1,3)}
          if (MA==2){MA=c(1,2)}else if (MA==3){MA=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        for (i in 1:nrow(models3)){#first step: 1-1,2-2,3
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(2,3)}
          if (AR==2){AR=c(1,2)}else if (AR==3){AR=c(2,3)}
          if (MA==2){MA=c(1,2)}else if (MA==3){MA=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        sep3.2=which(models3[,1]==2 | models3[,2]==2 | models3[,3]==2 | models3[,4]==2)
        if(length(sep3.2)!=0){
        models3.2=models3[sep3.2,] 
        
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,3
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,3)}
          if (ma==3){ma=c(1,3)}
          if (AR==3){AR=c(1,3)}
          if (MA==3){MA=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,2,3
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(1,2,3)}
          if (ma==3){ma=c(1,2,3)}
          if (AR==3){AR=c(1,2,3)}
          if (MA==3){MA=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        for (i in 1:nrow(models3.2)){#first step: 1-2-2,3
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if (ar==3){ar=c(2,3)}
          if (ma==3){ma=c(2,3)}
          if (AR==3){AR=c(2,3)}
          if (MA==3){MA=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        }
        l2 <- rep(list(0:2), n)
        possible2<-expand.grid(l2)
        for (i in 1:nrow(possible2)){
          order<-possible2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,2)}
            if(ma==1){ma=c(1,2)}
            if(AR==1){AR=c(1,2)}
            if(MA==1){MA=c(1,2)}
            mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        l3 <- rep(list(0:4), n)
        possible3<-expand.grid(l3)
        for (i in 1:nrow(possible3)){
          order<-possible3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,3)}else if(ar==2){ar=c(2,3)} else if(ar==4){ar=c(1,2,3)}
            if(ma==1){ma=c(1,3)}else if(ma==2){ma=c(2,3)} else if(ma==4){ma=c(1,2,3)}
            if(AR==1){AR=c(1,3)}else if(AR==2){AR=c(2,3)} else if(AR==4){AR=c(1,2,3)}
            if(MA==1){MA=c(1,3)}else if(MA==2){MA=c(2,3)} else if(MA==4){MA=c(1,2,3)}
            mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        print("first model up to order 3");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second model up to order 3")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third model up to order 3")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth model up to order 3")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth model up to order 3")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth model up to order 3")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh model up to order 3")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth model up to order 3")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth model up to order 3")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth model up to order 3")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
      }
      
      if (max.order==2){
        l <- rep(list(0:2), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]
         for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        sep2<-which(possible[,1]==2 | possible[,2]==2 | possible[,3]==2 | possible[,4]==2)
        models2<-possible[sep2,]
       
        for (i in 1:nrow(models2)){#first step 2
          order<-models2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          if (AR==2){AR=c(1,2)}
          if (MA==2){MA=c(1,2)}
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }}
        
        l2 <- rep(list(0:2), n)
        possible2<-expand.grid(l2)
        for (i in 1:nrow(possible2)){
          order<-possible2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,2)}
            if(ma==1){ma=c(1,2)}
            if(AR==1){AR=c(1,2)}
            if(MA==1){MA=c(1,2)}
            mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        
        print("first model up to order 2");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second model order 2")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third model order 2")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth model up to order 2")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth model up to order 2")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth model up to order 2")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh model up to order 2")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth model up to order 2")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth model up to order 2")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth model up to order 2")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
        
      }
      
      if (max.order==1){
        l <- rep(list(0:1), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]
         for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=unlist(order[3])
          MA=unlist(order[4])
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }
          }
          }
        }
        print("first model up to order 1");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        print("AR="); print(first$AR)
        print("MA="); print(first$MA)
        
        print("second model order 1")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        print("AR="); print(second$AR)
        print("MA="); print(second$MA)
        
        print("third model order 1")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        print("AR="); print(third$AR)
        print("MA="); print(third$MA)
        
        print("fourth model up to order 1")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        print("AR="); print(fourth$AR)
        print("MA="); print(fourth$MA)
        
        print("fifth model up to order 1")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        print("AR="); print(fifth$AR)
        print("MA="); print(fifth$MA)
        
        print("sixth model up to order 1")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        print("AR="); print(sixth$AR)
        print("MA="); print(sixth$MA)
        
        print("seventh model up to order 1")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        print("AR="); print(seventh$AR)
        print("MA="); print(seventh$MA)
        
        print("eighth model up to order 1")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        print("AR="); print(eighth$AR)
        print("MA="); print(eighth$MA)
        
        print("ninth model up to order 1")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        print("AR="); print(ninth$AR)
        print("MA="); print(ninth$MA)
        
        print("tenth model up to order 1")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        print("AR="); print(tenth$AR)
        print("MA="); print(tenth$MA)
      }
      #end seasonal models
    }else{#start non-seasonal models
      print(paste("Friedman Rank test. P-value=",friedman$Pval))
      
      print(paste("Excluding seasonal models on search for the best model."))
      n <- 2
      
      if (max.order==4){
        l <- rep(list(0:4), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]  
         
        for (i in 1:nrow(models1)){#first step: 1-2-3-4
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep2<-which(possible[,1]==2 | possible[,2]==2)
        models2<-possible[sep2,]
        
        for (i in 1:nrow(models2)){#first step: 1-1,2-3-4
          order<-models2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        sep2.2<-which(possible[,1]==2 & possible[,2]==2)
        models2.2<-possible[sep2.2,]
         
        for (i in 1:nrow(models2.2)){
          order<-models2.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          if (ar==2){ar=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        for (i in 1:nrow(models2.2)){
          order<-models2.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          if (ma==2){ma=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep3<-which(possible[,1]==3 | possible[,2]==3)
        models3<-possible[sep3,]
       
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,2,3-4
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,3-4
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        for (i in 1:nrow(models3)){#first step: 1-1,2-2,3-4
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        sep3.2=which(models3[,1]==2 | models3[,2]==2)
        if(length(sep3.2)!=0){
        models3.2=models3[sep3.2,]
        
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,3-4
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)}
          if (ma==3){ma=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,2,3-4
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)}
          if (ma==3){ma=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        for (i in 1:nrow(models3.2)){#first step: 1-2-2,3-4
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)}
          if (ma==3){ma=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        }
       
        sep4<-which(possible[,1]==4 | possible[,2]==4)
        models4<-possible[sep4,]
        
        #first step 1-2-3-1,4
        for (i in 1:nrow(models4)){#first step 1-2-3-1,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==4){ar=c(1,4)}
          if (ma==4){ma=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-3-2,4
        for (i in 1:nrow(models4)){#first step 1-2-3-2,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==4){ar=c(2,4)}
          if (ma==4){ma=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-3-3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==4){ar=c(3,4)}
          if (ma==4){ma=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-3-1,2,4
        for (i in 1:nrow(models4)){#first step 1-2-3-1,2,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==4){ar=c(1,2,4)}
          if (ma==4){ma=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-3-1,3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-1,3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==4){ar=c(1,3,4)}
          if (ma==4){ma=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-3-1,2,3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-1,2,3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==4){ar=c(1,2,3,4)}
          if (ma==4){ma=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-3-2,3,4
        for (i in 1:nrow(models4)){#first step 1-2-3-2,3,4
          order<-models4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==4){ar=c(2,3,4)}
          if (ma==4){ma=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep4.2=which(models4[,1]==2 | models4[,2]==2)
        if(length(sep4.2)!=0){
        models4.2=models4[sep4.2,]
        
        #first step 1-1,2-3-1,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-1,2-3-2,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-2,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-1,2-3-3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-1,2-3-1,2,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,2,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-1,2-3-1,3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-1,2-3-1,2,3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-1,2,3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-1,2-3-2,3,4
        for (i in 1:nrow(models4.2)){#first step 1-1,2-3-2,3,4
          order<-models4.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==4){ma=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        }
        sep4.3=which(models4[,1]==3 | models4[,2]==3)
        if(length(sep4.3)!=0){
        models4.3=models4[sep4.3,]
        
        #first step 1-2-1,2,3-1,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,2,3-2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,2,3-3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,2,3-1,2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,2,3-1,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,2,3-1,2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-1,2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,2,3-2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,2,3-2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,3-1,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,3-2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,3-3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,3-1,2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,3-1,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,3-1,2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-1,2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-1,3-2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-1,3-2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-2,3-1,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-2,3-2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-2,3-3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-2,3-1,2,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,2,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-2,3-1,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-2,3-1,2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-1,2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        #first step 1-2-2,3-2,3,4
        for (i in 1:nrow(models4.3)){#first step 1-2-2,3-2,3,4
          order<-models4.3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,3,4)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        sep4.3.2=which(models4.3[,1]==2 | models4.3[,2]==2)
        if(length(sep4.3.2)!=0){
        models4.3.2=models4.3[sep4.3.2,]#zero models
        
        #first step 1-1,2-1,2,3-1,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-1,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,2,3-2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        # #first step 1-1,2-1,2,3-3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,2,3-1,2,4
        for (i in 1:nrow(models4.3.2)){#first step  1-1,2-1,2,3-1,2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,2,3-1,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-1,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,2,3-1,2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-1,2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(1,2,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,2,3-2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,2,3-2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,2,3)} else if (ma==4){ma=c(2,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,3-1,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,3-2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,3-3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,3-1,2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,3-1,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,3-1,2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-1,2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(1,2,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-1,3-2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-1,3-2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(1,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(1,3)} else if (ma==4){ma=c(2,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-2,3-1,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-2,3-2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-2,3-3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-2,3-1,2,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,2,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-2,3-1,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-2,3-1,2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-1,2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(1,2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(1,2,3,4)}

          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        #first step 1-1,2-2,3-2,3,4
        for (i in 1:nrow(models4.3.2)){#first step 1-1,2-2,3-2,3,4
          order<-models4.3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)

          if (ar==2){ar=c(1,2)} else if (ar==3){ar=c(2,3)} else if (ar==4){ar=c(2,3,4)}
          if (ma==2){ma=c(1,2)} else if (ma==3){ma=c(2,3)} else if (ma==4){ma=c(2,3,4)}


          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          #print(mod$RMC)
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              third=second
              second=first
              first<-mod    ;       print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<second$mbic){
              third=second
              second<-mod;                  print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
            }else if  (mod$mbic<third$mbic){
              third<-mod;print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)}
          }
        }}
        }
        }
        
        l2 <- rep(list(0:2), n)
        possible2<-expand.grid(l2)
        for (i in 1:nrow(possible2)){
          order<-possible2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=0
          MA=0
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,2)}
            if(ma==1){ma=c(1,2)}
            mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        l3 <- rep(list(0:4), n)
        possible3<-expand.grid(l3)
        for (i in 1:nrow(possible3)){
          order<-possible3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=0
          MA=0
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,3)}else if(ar==2){ar=c(2,3)} else if(ar==4){ar=c(1,2,3)}
            if(ma==1){ma=c(1,3)}else if(ma==2){ma=c(2,3)} else if(ma==4){ma=c(1,2,3)}
           mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        l4 <- rep(list(0:8), n)
        possible4<-expand.grid(l4)
        for (i in 1:nrow(possible4)){
          order<-possible4[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=0
          MA=0
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,4)}else if(ar==2){ar=c(2,4)}else if(ar==3){ar=c(3,4)}else if(ar==5){ar=c(1,2,4)}else if(ar==6){ar=c(1,3,4)}else if(ar==7){ar=c(2,3,4)}else if(ar==8){ar=c(1,2,3,4)}
            if(ma==1){ma=c(1,4)}else if(ma==2){ma=c(2,4)}else if(ma==3){ma=c(3,4)}else if(ma==5){ma=c(1,2,4)}else if(ma==6){ma=c(1,3,4)}else if(ma==7){ma=c(2,3,4)}else if(ma==8){ma=c(1,2,3,4)}
           mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA=");print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        print("first model up to order 4");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        print("second model up to order 4")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        print("third model up to order 4")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        print("fourth model up to order 4")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        print("fifth model up to order 4")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        print("sixth model up to order 4")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        print("seventh model up to order 4")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        print("eighth model up to order 4")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        print("ninth model up to order 4")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        print("tenth model up to order 4")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
      }
      
      if (max.order==3){
        
        l <- rep(list(0:3), n)
        possible<-expand.grid(l)#15 models
        models1<-possible[3:nrow(possible),]
       
        for (i in 1:nrow(models1)){#first step: 1-2-3
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep2<-which(possible[,1]==2 | possible[,2]==2)
        models2<-possible[sep2,]
     
        for (i in 1:nrow(models2)){#first step: 1-1,2-3
          order<-models2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        sep2.2<-which(possible[,1]==2 & possible[,2]==2)
        models2.2<-possible[sep2.2,]
        
        for (i in 1:nrow(models2.2)){
          order<-models2.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          if (ar==2){ar=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        for (i in 1:nrow(models2.2)){
          order<-models2.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          if (ma==2){ma=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep3<-which(possible[,1]==3 | possible[,2]==3)
        models3<-possible[sep3,]
        
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,2,3
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        for (i in 1:nrow(models3)){#first step: 1-1,2-1,3
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(1,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        for (i in 1:nrow(models3)){#first step: 1-1,2-2,3
          order<-models3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==2){ar=c(1,2)}else if (ar==3){ar=c(2,3)}
          if (ma==2){ma=c(1,2)}else if (ma==3){ma=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep3.2=which(models3[,1]==2 | models3[,2]==2)
       
        if(length(sep3.2)!=0){
        models3.2=models3[sep3.2,] 
       
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,3
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,3)}
          if (ma==3){ma=c(1,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        for (i in 1:nrow(models3.2)){#first step: 1-2-1,2,3
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(1,2,3)}
          if (ma==3){ma=c(1,2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
       
        
        
        for (i in 1:nrow(models3.2)){#first step: 1-2-2,3
          order<-models3.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          if (ar==3){ar=c(2,3)}
          if (ma==3){ma=c(2,3)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        }
        
        l2 <- rep(list(0:2), n)
        possible2<-expand.grid(l2)
        for (i in 1:nrow(possible2)){
          order<-possible2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=0
          MA=0
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,2)}
            if(ma==1){ma=c(1,2)}
            mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        l3 <- rep(list(0:4), n)
        possible3<-expand.grid(l3)
        for (i in 1:nrow(possible3)){
          order<-possible3[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=0
          MA=0
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,3)}else if(ar==2){ar=c(2,3)} else if(ar==4){ar=c(1,2,3)}
            if(ma==1){ma=c(1,3)}else if(ma==2){ma=c(2,3)} else if(ma==4){ma=c(1,2,3)}
            mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        print("first model up to order 3");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        print("second model up to order 3")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        print("third model up to order 3")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        print("fourth model up to order 3")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        print("fifth model up to order 3")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        print("sixth model up to order 3")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        print("seventh model up to order 3")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        print("eighth model up to order 3")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        print("ninth model up to order 3")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        print("tenth model up to order 3")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
      }
      
      if (max.order==2){
        l <- rep(list(0:2), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]#7 models
        
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep2<-which(possible[,1]==2 | possible[,2]==2)
        models2<-possible[sep2,]#5 models to set c(1,2) order
        
        for (i in 1:nrow(models2)){#first step 2
          order<-models2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          if (ar==2){ar=c(1,2)}
          if (ma==2){ma=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        sep2.2<-which(possible[,1]==2 & possible[,2]==2)
        models2.2<-possible[sep2.2,]
          
        for (i in 1:nrow(models2.2)){
          order<-models2.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          if (ar==2){ar=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        for (i in 1:nrow(models2.2)){
          order<-models2.2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          if (ma==2){ma=c(1,2)}
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        l2 <- rep(list(0:2), n)
        possible2<-expand.grid(l2)
        for (i in 1:nrow(possible2)){
          order<-possible2[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=0
          MA=0
          
          if(length(which(c(ar,ma,AR,MA)==0))<=2){
            if(ar==1){ar=c(1,2)}
            if(ma==1){ma=c(1,2)}
            mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
            if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
              if (mod$mbic<first$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second=first;
                first<-mod;                      
                print("new first"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<second$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third=second;
                second<-mod;                  
                print("new second"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<third$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth=third;
                third<-mod;
                print("new third"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fourth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth=fourth;
                fourth<-mod;
                print("new fourth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<fifth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth=fifth;
                fifth<-mod;
                print("new fifth");print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<sixth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh=sixth;
                sixth<-mod;
                print("new sixth"); print("ar="); print(ar); print("ma="); print(ma); print("AR=");print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<seventh$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth=seventh;
                seventh<-mod;
                print("new seventh"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<eighth$mbic){
                tenth=ninth;
                ninth=eighth;
                eighth<-mod;
                print("new eighth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<ninth$mbic){
                tenth=ninth;
                ninth<-mod;
                print("new ninth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }else if  (mod$mbic<tenth$mbic){
                tenth<-mod;
                print("new tenth"); print("ar="); print(ar); print("ma="); print(ma); print("AR="); print(AR); print("MA="); print(MA)
              }
            }
            }
          }
        }
        
        
        print("first model up to order 2");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        print("second model up to order 2")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        print("third model up to order 2")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        print("fourth model up to order 2")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        print("fifth model up to order 2")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        print("sixth model up to order 2")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        print("seventh model up to order 2")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        print("eighth model up to order 2")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        print("ninth model up to order 2")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        print("tenth model up to order 2")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
      }
      
      if (max.order==1){
        l <- rep(list(0:1), n)
        possible<-expand.grid(l)
        models1<-possible[3:nrow(possible),]#2 models
         
        for (i in 1:nrow(models1)){#first level 1
          order<-models1[i,]
          ar=unlist(order[1])
          ma=unlist(order[2])
          AR=c(0)
          MA=c(0)
          
          mod<-try( mkarma(y=y,ar=ar,ma=ma,AR=AR,MA=MA,exvar=exvar,steps=steps,graph=graph,print=print,S=S,resid=resid))
          if (mod$RMC==0 ){if(all(mod$pvalues<=0.05)& mod$p_ljungbox>=0.05 & mod$p_jarquebera>=0.05 & mod$p_arch>=0.05 &all(mod$roots>=1)){
            if (mod$mbic<first$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second=first;
              first<-mod;                      
              print("new first"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<second$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third=second;
              second<-mod;                  
              print("new second"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<third$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth=third;
              third<-mod;
              print("new third"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fourth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth=fourth;
              fourth<-mod;
              print("new fourth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<fifth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth=fifth;
              fifth<-mod;
              print("new fifth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<sixth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh=sixth;
              sixth<-mod;
              print("new sixth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<seventh$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth=seventh;
              seventh<-mod;
              print("new seventh"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<eighth$mbic){
              tenth=ninth;
              ninth=eighth;
              eighth<-mod;
              print("new eighth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<ninth$mbic){
              tenth=ninth;
              ninth<-mod;
              print("new ninth"); print("ar="); print(ar); print("ma="); print(ma)
            }else if  (mod$mbic<tenth$mbic){
              tenth<-mod;
              print("new tenth"); print("ar="); print(ar); print("ma="); print(ma)
            }
          }
          }}
        
        print("first model up to order 1");print(first$mbic)
        print("ar="); print(first$ar)
        print("ma="); print(first$ma)
        
        print("second model up to order 1")
        print("ar="); print(second$ar)
        print("ma="); print(second$ma)
        
        print("third model up to order 1")
        print("ar="); print(third$ar)
        print("ma="); print(third$ma)
        
        print("fourth model up to order 1")
        print("ar="); print(fourth$ar)
        print("ma="); print(fourth$ma)
        
        print("fifth model up to order 1")
        print("ar="); print(fifth$ar)
        print("ma="); print(fifth$ma)
        
        print("sixth model up to order 1")
        print("ar="); print(sixth$ar)
        print("ma="); print(sixth$ma)
        
        print("seventh model up to order 1")
        print("ar="); print(seventh$ar)
        print("ma="); print(seventh$ma)
        
        print("eighth model up to order 1")
        print("ar="); print(eighth$ar)
        print("ma="); print(eighth$ma)
        
        print("ninth model up to order 1")
        print("ar="); print(ninth$ar)
        print("ma="); print(ninth$ma)
        
        print("tenth model up to order 1")
        print("ar="); print(tenth$ar)
        print("ma="); print(tenth$ma)
        
      }
    }#non-seasonal models
  }
  
  first$sortmbic="first"
  second$sortmbic="second"
  third$sortmbic="third"
  fourth$sortmbic="fourth"
  fifth$sortmbic="fifth"
  sixth$sortmbic="sixth"
  seventh$sortmbic="seventh"
  eighth$sortmbic="eighth"
  ninth$sortmbic="ninth"
  tenth$sortmbic="tenth"
  
  all=c(first$accuracyfitted[,4],second$accuracyfitted[,4], third$accuracyfitted[,4], fourth$accuracyfitted[,4],fifth$accuracyfitted[,4],sixth$accuracyfitted[,4],seventh$accuracyfitted[,4],eighth$accuracyfitted[,4],ninth$accuracyfitted[,4],tenth$accuracyfitted[,4])
  
  compare=sort(all)
  if(first$accuracyfitted[,4]==compare[1]){m.first=first; m.first$sortmape="first"}else
    if(second$accuracyfitted[,4]==compare[1]){m.first=second; m.first$sortmape="first"}else
      if(third$accuracyfitted[,4]==compare[1]){m.first=third; m.first$sortmape="first"}else
        if(fourth$accuracyfitted[,4]==compare[1]){m.first=fourth; m.first$sortmape="first"}else
          if(fifth$accuracyfitted[,4]==compare[1]){m.first=fifth; m.first$sortmape="first"}else
            if(sixth$accuracyfitted[,4]==compare[1]){m.first=sixth; m.first$sortmape="first"}else
              if(seventh$accuracyfitted[,4]==compare[1]){m.first=seventh; m.first$sortmape="first"}else
                if(eighth$accuracyfitted[,4]==compare[1]){m.first=eighth; m.first$sortmape="first"}else
                  if(ninth$accuracyfitted[,4]==compare[1]){m.first=ninth; m.first$sortmape="first"}else
                    if(tenth$accuracyfitted[,4]==compare[1]){m.first=tenth; m.first$sortmape="first"}
  if(first$accuracyfitted[,4]==compare[2]){m.second=first; m.second$sortmape="second"}else
    if(second$accuracyfitted[,4]==compare[2]){m.second=second; m.second$sortmape="second"}else
      if(third$accuracyfitted[,4]==compare[2]){m.second=third; m.second$sortmape="second"}else
        if(fourth$accuracyfitted[,4]==compare[2]){m.second=fourth; m.second$sortmape="second"}else
          if(fifth$accuracyfitted[,4]==compare[2]){m.second=fifth; m.second$sortmape="second"}else
            if(sixth$accuracyfitted[,4]==compare[2]){m.second=sixth; m.second$sortmape="second"}else
              if(seventh$accuracyfitted[,4]==compare[2]){m.second=seventh; m.second$sortmape="second"}else
                if(eighth$accuracyfitted[,4]==compare[2]){m.second=eighth; m.second$sortmape="second"}else
                  if(ninth$accuracyfitted[,4]==compare[2]){m.second=ninth; m.second$sortmape="second"}else
                    if(tenth$accuracyfitted[,4]==compare[2]){m.second=tenth; m.second$sortmape="second"}
  if(first$accuracyfitted[,4]==compare[3]){m.third=first; m.third$sortmape="third"}else
    if(second$accuracyfitted[,4]==compare[3]){m.third=second; m.third$sortmape="third"}else
      if(third$accuracyfitted[,4]==compare[3]){m.third=third; m.third$sortmape="third"}else
        if(fourth$accuracyfitted[,4]==compare[3]){m.third=fourth; m.third$sortmape="third"}else
          if(fifth$accuracyfitted[,4]==compare[3]){m.third=fifth; m.third$sortmape="third"}else
            if(sixth$accuracyfitted[,4]==compare[3]){m.third=sixth; m.third$sortmape="third"}else
              if(seventh$accuracyfitted[,4]==compare[3]){m.third=seventh; m.third$sortmape="third"}else
                if(eighth$accuracyfitted[,4]==compare[3]){m.third=eighth; m.third$sortmape="third"}else
                  if(ninth$accuracyfitted[,4]==compare[3]){m.third=ninth; m.third$sortmape="third"}else
                    if(tenth$accuracyfitted[,4]==compare[3]){m.third=tenth; m.third$sortmape="third"}
  
  if(first$accuracyfitted[,4]==compare[4]){m.fourth=first; m.fourth$sortmape="fourth"}else
    if(second$accuracyfitted[,4]==compare[4]){m.fourth=second; m.fourth$sortmape="fourth"}else
      if(third$accuracyfitted[,4]==compare[4]){m.fourth=third; m.fourth$sortmape="fourth"}else
        if(fourth$accuracyfitted[,4]==compare[4]){m.fourth=fourth; m.fourth$sortmape="fourth"}else
          if(fifth$accuracyfitted[,4]==compare[4]){m.fourth=fifth; m.fourth$sortmape="fourth"}else
            if(sixth$accuracyfitted[,4]==compare[4]){m.fourth=sixth; m.fourth$sortmape="fourth"}else
              if(seventh$accuracyfitted[,4]==compare[4]){m.fourth=seventh; m.fourth$sortmape="fourth"}else
                if(eighth$accuracyfitted[,4]==compare[4]){m.fourth=eighth; m.fourth$sortmape="fourth"}else
                  if(ninth$accuracyfitted[,4]==compare[4]){m.fourth=ninth; m.fourth$sortmape="fourth"}else
                    if(tenth$accuracyfitted[,4]==compare[4]){m.fourth=tenth; m.fourth$sortmape="fourth"}
  
  if(first$accuracyfitted[,4]==compare[5]){m.fifth=first; m.fifth$sortmape="fifth"}else
    if(second$accuracyfitted[,4]==compare[5]){m.fifth=second; m.fifth$sortmape="fifth"}else
      if(third$accuracyfitted[,4]==compare[5]){m.fifth=third; m.fifth$sortmape="fifth"}else
        if(fourth$accuracyfitted[,4]==compare[5]){m.fifth=fourth; m.fifth$sortmape="fifth"}else
          if(fifth$accuracyfitted[,4]==compare[5]){m.fifth=fifth; m.fifth$sortmape="fifth"}else
            if(sixth$accuracyfitted[,4]==compare[5]){m.fifth=sixth; m.fifth$sortmape="fifth"}else
              if(seventh$accuracyfitted[,4]==compare[5]){m.fifth=seventh; m.fifth$sortmape="fifth"}else
                if(eighth$accuracyfitted[,4]==compare[5]){m.fifth=eighth; m.fifth$sortmape="fifth"}else
                  if(ninth$accuracyfitted[,4]==compare[5]){m.fifth=ninth; m.fifth$sortmape="fifth"}else
                    if(tenth$accuracyfitted[,4]==compare[5]){m.fifth=tenth; m.fifth$sortmape="fifth"}
  
  if(first$accuracyfitted[,4]==compare[6]){m.sixth=first; m.sixth$sortmape="sixth"}else
    if(second$accuracyfitted[,4]==compare[6]){m.sixth=second; m.sixth$sortmape="sixth"}else
      if(third$accuracyfitted[,4]==compare[6]){m.sixth=third; m.sixth$sortmape="sixth"}else
        if(fourth$accuracyfitted[,4]==compare[6]){m.sixth=fourth; m.sixth$sortmape="sixth"}else
          if(fifth$accuracyfitted[,4]==compare[6]){m.sixth=fifth; m.sixth$sortmape="sixth"}else
            if(sixth$accuracyfitted[,4]==compare[6]){m.sixth=sixth; m.sixth$sortmape="sixth"}else
              if(seventh$accuracyfitted[,4]==compare[6]){m.sixth=seventh; m.sixth$sortmape="sixth"}else
                if(eighth$accuracyfitted[,4]==compare[6]){m.sixth=eighth; m.sixth$sortmape="sixth"}else
                  if(ninth$accuracyfitted[,4]==compare[6]){m.sixth=ninth; m.sixth$sortmape="sixth"}else
                    if(tenth$accuracyfitted[,4]==compare[6]){m.sixth=tenth; m.sixth$sortmape="sixth"}
  
  if(first$accuracyfitted[,4]==compare[7]){m.seventh=first; m.seventh$sortmape="seventh"}else
    if(second$accuracyfitted[,4]==compare[7]){m.seventh=second; m.seventh$sortmape="seventh"}else
      if(third$accuracyfitted[,4]==compare[7]){m.seventh=third; m.seventh$sortmape="seventh"}else
        if(fourth$accuracyfitted[,4]==compare[7]){m.seventh=fourth; m.seventh$sortmape="seventh"}else
          if(fifth$accuracyfitted[,4]==compare[7]){m.seventh=fifth; m.seventh$sortmape="seventh"}else
            if(sixth$accuracyfitted[,4]==compare[7]){m.seventh=sixth; m.seventh$sortmape="seventh"}else
              if(seventh$accuracyfitted[,4]==compare[7]){m.seventh=seventh; m.seventh$sortmape="seventh"}else
                if(eighth$accuracyfitted[,4]==compare[7]){m.seventh=eighth; m.seventh$sortmape="seventh"}else
                  if(ninth$accuracyfitted[,4]==compare[7]){m.seventh=ninth; m.seventh$sortmape="seventh"}else
                    if(tenth$accuracyfitted[,4]==compare[7]){m.seventh=tenth; m.seventh$sortmape="seventh"}
  
  if(first$accuracyfitted[,4]==compare[8]){m.eighth=first; m.eighth$sortmape="eighth"}else
    if(second$accuracyfitted[,4]==compare[8]){m.eighth=second; m.eighth$sortmape="eighth"}else
      if(third$accuracyfitted[,4]==compare[8]){m.eighth=third; m.eighth$sortmape="eighth"}else
        if(fourth$accuracyfitted[,4]==compare[8]){m.eighth=fourth; m.eighth$sortmape="eighth"}else
          if(fifth$accuracyfitted[,4]==compare[8]){m.eighth=fifth; m.eighth$sortmape="eighth"}else
            if(sixth$accuracyfitted[,4]==compare[8]){m.eighth=sixth; m.eighth$sortmape="eighth"}else
              if(seventh$accuracyfitted[,4]==compare[8]){m.eighth=seventh; m.eighth$sortmape="eighth"}else
                if(eighth$accuracyfitted[,4]==compare[8]){m.eighth=eighth; m.eighth$sortmape="eighth"}else
                  if(ninth$accuracyfitted[,4]==compare[8]){m.eighth=ninth; m.eighth$sortmape="eighth"}else                    
                    if(tenth$accuracyfitted[,4]==compare[8]){m.eighth=tenth; m.eighth$sortmape="eighth"}
  
  if(first$accuracyfitted[,4]==compare[9]){m.ninth=first; m.ninth$sortmape="ninth"}else
    if(second$accuracyfitted[,4]==compare[9]){m.ninth=second; m.ninth$sortmape="ninth"}else
      if(third$accuracyfitted[,4]==compare[9]){m.ninth=third; m.ninth$sortmape="ninth"}else
        if(fourth$accuracyfitted[,4]==compare[9]){m.ninth=fourth; m.ninth$sortmape="ninth"}else
          if(fifth$accuracyfitted[,4]==compare[9]){m.ninth=fifth; m.ninth$sortmape="ninth"}else
            if(sixth$accuracyfitted[,4]==compare[9]){m.ninth=sixth; m.ninth$sortmape="ninth"}else
              if(seventh$accuracyfitted[,4]==compare[9]){m.ninth=seventh; m.ninth$sortmape="ninth"}else
                if(eighth$accuracyfitted[,4]==compare[9]){m.ninth=eighth; m.ninth$sortmape="ninth"}else
                  if(ninth$accuracyfitted[,4]==compare[9]){m.ninth=ninth; m.ninth$sortmape="ninth"}else
                    if(tenth$accuracyfitted[,4]==compare[9]){m.ninth=tenth; m.ninth$sortmape="ninth"}
  
  if(first$accuracyfitted[,4]==compare[10]){m.tenth=first; m.tenth$sortmape="tenth"}else   
    if(second$accuracyfitted[,4]==compare[10]){m.tenth=second; m.tenth$sortmape="tenth"}else
      if(third$accuracyfitted[,4]==compare[10]){m.tenth=third; m.tenth$sortmape="tenth"}else
        if(fourth$accuracyfitted[,4]==compare[10]){m.tenth=fourth; m.tenth$sortmape="tenth"}else
          if(fifth$accuracyfitted[,4]==compare[10]){m.tenth=fifth; m.tenth$sortmape="tenth"}else
            if(sixth$accuracyfitted[,4]==compare[10]){m.tenth=sixth; m.tenth$sortmape="tenth"}else
              if(seventh$accuracyfitted[,4]==compare[10]){m.tenth=seventh; m.tenth$sortmape="tenth"}else
                if(eighth$accuracyfitted[,4]==compare[10]){m.tenth=eighth; m.tenth$sortmape="tenth"}else
                  if(ninth$accuracyfitted[,4]==compare[10]){m.tenth=ninth; m.tenth$sortmape="tenth"}else
                    if(tenth$accuracyfitted[,4]==compare[10]){m.tenth=tenth; m.tenth$sortmape="tenth"}
  toc()
  sink()
  if(sort=="MAPE"){ print("sort=MAPE so the printed 10 lowest mbics were sorted by MAPE on the return")
    return(list(first=m.first,second=m.second,third=m.third,fourth=m.fourth,fifth=m.fifth,sixth=m.sixth,seventh=m.seventh,eighth=m.eighth,ninth=m.ninth,tenth=m.tenth))
  }else{print("sort=mbic 10 lowest mbics")
    return(list(first=first,second=second,third=third,fourth=fourth,fifth=fifth,sixth=sixth,seventh=seventh,eighth=eighth,ninth=ninth,tenth=tenth))
  }
  
}
