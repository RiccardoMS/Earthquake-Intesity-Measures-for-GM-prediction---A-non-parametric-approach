#working directory
setwd("C:/Users/aughi/Desktop/Progetto Stat/Dataset_Southern_Italy")

#import dataset
library(readr)
DATA_0.50_NO <- read.csv(" FAS_000_50_NoOutlier.csv")
head(DATA_0.50_NO)
dim(DATA_0.50_NO)
colnames(DATA_0.50_NO)
DATA_0.50_NO<-DATA_0.50_NO[,-1]
attach(DATA_0.50_NO)

#to define the model, consider only the values for which Usable VS30 is defined
detach(DATA_0.50_NO)
DATA_0.50_NO<-DATA_0.50_NO[-which(is.na(DATA_0.50_NO$UsableVS30)),]
attach(DATA_0.50_NO)

#ML-logY (Dipo bins)
range(Dipo)
Dipo.lab=c("0-50 [km]","50-100 [km]","100-150 [km]")
Dipo.bin <- seq(min(Dipo)-0.0001, max(Dipo)+0.0001, length.out = 4)
Dipo.seq<-cut(Dipo,Dipo.bin)
Dipo.col<- rep(NA, length(Dipo))
for (i in 1:length(Dipo)){
  if((Dipo.seq[i]==levels(Dipo.seq)[1])==T)
    Dipo.col[i]='red'
  if((Dipo.seq[i]==levels(Dipo.seq)[2])==T)
    Dipo.col[i]='green'
  if((Dipo.seq[i]==levels(Dipo.seq)[3])==T)
    Dipo.col[i]='blue'
}
plot(ML,LogY, main='Response Against ML',xlab='Magnitude',ylab='Log10 of GeomMean',col=Dipo.col)
legend('bottomright',Dipo.lab,col=c('red','green','blue'),pch='o')

#regressione lineare a tratti?
library(segmented)
mod<-lm(LogY~ML)
piecewise.lin.mod<-segmented(mod, psi=5) #modFM
plot(ML,LogY,col='green')
points(ML,piecewise.lin.mod$fitted.values,cex=0.8)
summary(piecewise.lin.mod)
plot(ML,piecewise.lin.mod$fitted.values)
segments(ML,piecewise.lin.mod$fitted.values,ML,LogY,col=Dipo.col,cex=0.1)
plot(ML,piecewise.lin.mod$residuals,col=Dipo.col,cex=0.1)
segments(ML,piecewise.lin.mod$residuals,ML,seq(0,0,length.out = length(ML)),col=Dipo.col)


#Dipo-LogY (ML bins)
range(ML)
ML.labels=c("2.5-3","3-3.5","3.5-4","4-4.5","4.5-5","5-5.5")
ML.bin <- seq(min(ML)-0.01, max(ML)+0.01, length.out = 7)
ML.seq<-cut(ML,ML.bin)
ML.col<- rep(NA, length(ML))
col.span<-rainbow(6)
for (i in 1:length(ML)){
  for (j in 1:6){
    if((ML.seq[i]==levels(ML.seq)[j])==T)
      ML.col[i]=col.span[j]
  }
}
plot(Dipo,LogY, main='Attenuation of Responde with Dipo',xlab='Dipo[km]',ylab='Log10 of GeomMean',col=ML.col)
legend('bottomright',ML.labels,col=col.span,pch='o',title='Magnitude Level')
#to evaluate:possible transofrmation of DIpo in log(DIpo)

library(minpack.lm)

#FUNZIONE FITTIZIA

func<-function(const,a1,a2,mh,c1,c3,mref,h){
 pred=const+ a1*ifelse((ML- mh)<=0,ML-mh,0)+ a2*ifelse((ML- mh)>0,ML-mh,0)+ c1*(ML-mref)*log10(sqrt(Depi^2+h^2))+ 0*log10(sqrt(Depi^2+h^2)) + c3*sqrt(Depi^2+h^2)
 par(mfrow=c(2,2))
 plot(ML,LogY)
 plot(ML,pred,cex=0.6,col='red')
 plot(Dipo,LogY)
 plot(Dipo,pred,cex=0.6,col='red')
}

func(-2,0.55,0.65,3.848,0.6,-0.003,3.8,3.3)


#SENZA C2
MODELLO1     <-nlsLM(LogY~const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
             + a2*ifelse((ML- mh)>0,ML-mh,0) 
             +c1*(ML-mref)*log10(sqrt(Depi^2+h^2))
             +0*log10(sqrt(Depi^2+h^2))
             +c3* sqrt(Depi^2+h^2) 
             +d1*log10(UsableVS30/800), 
             start=list(const=-1,a1=0.55, a2=0.65, mh=3.8, mref=4, h=7,c1=0.6,c3=-0.003,d1=-0.3), control = nls.lm.control(maxiter = 1000))
MODELLO1
### STABLE, BUT PAREMETERS MAY BE OVERESTIMATED

#RUMOR ON DEPI, terms with c1-c3
x=rnorm(n=1000,m=0,sd=1)
DepiF=Depi+sample(x,length(Depi),replace = T)
MODELLO2<-nlsLM(LogY~const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
                      + a2*ifelse((ML- mh)>0,ML-mh,0) 
                      +c1*(ML-mref)*log10(sqrt(DepiF^2+h^2))
                      +c2*log10(sqrt(Depi^2+h^2))
                      +c3* sqrt(DepiF^2+h^2) 
                      +d1*log10(UsableVS30/800), 
                      start=list(const=-1,a1=0.55, a2=0.65, mh=3.8, mref=4, h=7,c1=0.6,c2=-0.01,c3=-0.003,d1=-0.3), control = nls.lm.control(maxiter = 1000))
MODELLO2
###UNSTABLE, DIFF RUMORS AFFECTS DEEPLY ESTIMATE OF MREF; DON'T KNOW THE ERROR I'M INTRODUCING



#RUMOR ON DEPI, terms with c2
x=rnorm(n=1000,m=0,sd=1)
DepiF=Depi+sample(x,length(Depi),replace = T)
MODELLO3<-nlsLM(LogY~const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
                + a2*ifelse((ML- mh)>0,ML-mh,0) 
                +c1*(ML-mref)*log10(sqrt(Depi^2+h^2))
                +c2*log10(sqrt(DepiF^2+h^2))
                +c3* sqrt(Depi^2+h^2) 
                +d1*log10(UsableVS30/800), 
                start=list(const=-1,a1=0.55, a2=0.65, mh=3.8, mref=4, h=7,c1=0.6,c2=-0.01,c3=-0.003,d1=-0.3), control = nls.lm.control(maxiter = 1000))
MODELLO3
###APPEARS TO BE MORE STABLE

#NO RUMOR, FIXED EFFECTS OF C2
MODELLO4<-nlsLM(LogY~const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
                + a2*ifelse((ML- mh)>0,ML-mh,0) 
                +c1*(ML-mref)*log10(sqrt(Depi^2+h^2))
                -1.8*log10(sqrt(Depi^2+h^2))
                +c3* sqrt(Depi^2+h^2) 
                +d1*log10(UsableVS30/800), 
                start=list(const=-1,a1=0.55, a2=0.65, mh=3.8, mref=4, h=7,c1=0.6,c3=-0.003,d1=-0.3), control = nls.lm.control(maxiter = 1000))
MODELLO4
###GOOD IDEA; DEFINE A RANGE OF VARIATIONS FOR FIXED C2 AND FIND THE BEST MODEL
variations=seq(from = -1,to=-3,length.out=2000)
REStoMIN=rep(NA,length(variations))
for(i in 1:length(variations)){
  FictitiousMod<-nlsLM(LogY~const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
                       + a2*ifelse((ML- mh)>0,ML-mh,0) 
                       +c1*(ML-mref)*log10(sqrt(Depi^2+h^2))
                       +variations[i]*log10(sqrt(Depi^2+h^2))
                       +c3* sqrt(Depi^2+h^2) 
                       +d1*log10(UsableVS30/800), 
                       start=list(const=1,a1=0.55, a2=0.65, mh=3.8, mref=4, h=7,c1=0.6,c3=-0.003,d1=-0.3), control = nls.lm.control(maxiter = 1000))
  REStoMIN[i]=sum(abs(FictitiousMod$m$resid()))
  remove(FictitiousMod)
}
fixedC2=variations[which(REStoMIN==min(REStoMIN))] #NOTE: ONCE DEFINED THE OPTIMAL C2, IT IS A POSSIBLE 
                                                   #      IDEA TO SHRINK THE INTERVAL (I.E. 100 values of C2
                                                   #      to evaluate between 1.79 and 1.82) TO BOOST PRECISION

#--------->NO RUMOR, FIXED EFFECTS OF C2, BEST C2 BETWEEN (-1;-3)
MODELLO5<-nlsLM(LogY~const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
                + a2*ifelse((ML- mh)>0,ML-mh,0) 
                +c1*(ML-mref)*log10(sqrt(Depi^2+h^2))
                +fixedC2*log10(sqrt(Depi^2+h^2))
                +c3* sqrt(Depi^2+h^2) 
                +d1*log10(UsableVS30/800), 
                start=list(const=-1,a1=0.55, a2=0.65, mh=3.8, mref=4, h=7,c1=0.6,c3=-0.003,d1=-0.3), control = nls.lm.control(maxiter = 1000))
MODELLO5
par(mfrow=c(1,2))
scatter3D(ML, Depi, LogY, pch = 18, cex = 0.25, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY")
scatter3D(ML, Depi, MODELLO5$m$fitted(), pch = 18, cex = 0.25, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY")
par(mfrow=c(1,1))
scatter3D(ML, Depi,abs( MODELLO5$m$fitted()-LogY), pch = 18, cex = 0.25, 
          theta = 270, phi = 15, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY")
segments3D(ML, Depi,abs(MODELLO5$m$fitted()-LogY),ML, Depi,seq(from=0,to=0,length.out=length(Depi)))


############################ SAME PROCEDURE TO REPEAT ON DIPO #########################################




