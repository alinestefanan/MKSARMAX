#GENERATE MKSARMA SAMPLE 
source("sample.mkarma.R")
y<-sample.mkarma(150, beta=0.5290,phi=0.4229,theta=0,PHI=0.2971,THETA=0,alpha=15.6063)

#MKSARMAX BEST MODEL
source("auto.mkarma.R")
MKB=auto.mkarma(y,steps=12,max.order=1,type="all")
MKB$first
MKB$second

#MKSARMA FIT
source("fit.R")
MK=mkarma(y,ar=c(1),ma=c(0), AR=c(1),MA=c(0),S=12,steps=0,graph=T,print=T)


#MKSARMAX FIT
# tendency
t <- 1:length(y)
# deterministic seasonality
C<-cos(2*pi*t/12)

MKX=mkarma(y,steps=12,ar=c(1),ma=c(1), AR=c(1),MA=c(1),S=12,exvar=C,print=T,graph=T)
