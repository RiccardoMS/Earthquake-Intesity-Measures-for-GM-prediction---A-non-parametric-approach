##set working directory
setwd("~/Desktop/Progetto/dataset_cleaned")

##import DATA
library(readr)
library(lmerTest)
DATA<- read.csv("FAS_000_50_NoOutlier.csv")
load("~/Desktop/Progetto/dataset_cleaned/coeff_estimates_000.50HZ.RData")

#verify DATA structure
head(DATA)

#pre processing
#remove missing values of VS30
DATA<-DATA[-which(is.na(DATA$UsableVS30)),]

#define groups
EV.fact<-as.factor(DATA$ID)
ST.fact<-as.factor(DATA$Stat)

#set nonlinear coefficients
mh     <-coeff.estimates.000_50HZ[[2]][1]
mref   <-coeff.estimates.000_50HZ[[2]][2]
h      <-coeff.estimates.000_50HZ[[2]][3]

#MODEL:independent  mixed effects due to station and event
# LogY =  const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
#        + a2*ifelse((ML- mh)>0,ML-mh,0) 
#        +(c1*(ML-mref)+c2)*log10(sqrt(Dipo^2+h^2))
#        +c3* sqrt(Dipo^2+h^2)
#        +d1*log10(UsableVS30/800)
#        +randomeffectdependingonEV
#        +randomeffectdependingonST
## NOTE: we are supposing of an effect directly on the response value,non in one of the regressors
#        --> (1|EV.fact) : practically we account for a "casual intercept" depending on single event
#load package
library(lme4)
library(nlme)
attach(DATA)

#define DESIGN MATRIX
Z<-data.frame(X0=rep(1,dim(DATA)[1]),X1=ifelse((ML- mh)<=0,ML-mh,0),X2=ifelse((ML- mh)>0,ML-mh,0),
         X3=(ML-mref)*log10(sqrt(Dipo^2+h^2)),X4=log10(sqrt(Dipo^2+h^2)),X5=sqrt(Dipo^2+h^2),
         X6=log10(UsableVS30/800))
attach(Z)

##simple model, aov, diagnostic
# #LME.fit<-lmer(LogY~ X1+X2+X3+X4+X5+X6+(1|EV.fact)+(1|ST.fact), 
#               data=cbind(Z,EV.fact,ST.fact,LogY))
LME.fit<-lmer(LogY~ X1+X2+X3+X4+X5+X6+(1|ST.fact), 
              data=cbind(Z,EV.fact,ST.fact,LogY))
summary(LME.fit)

plot(LME.fit,xlab='fitted',ylab='residuals') #omoschedastici
qqnorm(residuals(LME.fit))
qqline(residuals(LME.fit),col='red') #gaussian

    # LME.fit1<-lme(LogY~ X1+X2+X3+X4+X5+X6, random=~1|EV.fact / ST.fact,data=cbind(Z,EV.fact,ST.fact,LogY))
  # summary(LME.fit1)
  # aov.lme<-aov(LME.fit1)
  # summary(aov.lme)

##Model with EC8
detach(Z)
Zdumm<-cbind.data.frame(Z,D1=DummyA,D2=DummyAStar,D3=DummyB,D4=DummyBStar,D5=DummyC,D6=DummyCStar,D7=DummyE)
attach(Zdumm)
LME.dumm<-lmer(LogY~X1+X2+X3+X4+X5+X6 
               +D1+D2+D3+D4+D5+D6+D7+(1|ST.fact), 
               data=cbind(Zdumm,EV.fact,ST.fact,LogY))


summary(LME.dumm) #poco significativi
vcov(LME.dumm)
plot(LME.dumm) #omoschedastici
          qqnorm(residuals(LME.dumm))
qqline(residuals(LME.dumm),col='red') #gaussian



#aov gives convergence error
# LME.dumm1<-lme(LogY~X1+X2+X3+X4+X5+X6+D1+D2+D3+D4+D5+D6+D7, random=~1|EV.fact / ST.fact,data=cbind(Zdumm,EV.fact,ST.fact,LogY))
# summary(LME.dumm1) 
# aov.lme.dumm<-aov(LME.dumm1)
# summary(aov.lme.dumm)

##
##X1:D1+X1:D2+X1:D3+X1:D4+X1:D5+X1:D6+X1:D7
##+X2:D1+X2:D2+X2:D3+X2:D4+X2:D5+X2:D6+X2:D7+
 # X3:D1+X3:D2+X3:D3+X3:D4+X3:D5+X3:D6+X3:D7+
#  X4:D1+X4:D2+X4:D3+X4:D4+X4:D5+X4:D6+X4:D7+
#  X5:D1+X5:D2+X5:D3+X5:D4+X5:D5+X5:D6+X5:D7
#+X6:D1+X6:D2+X6:D3+X6:D4+X6:D5+X6:D6+X6:D7
