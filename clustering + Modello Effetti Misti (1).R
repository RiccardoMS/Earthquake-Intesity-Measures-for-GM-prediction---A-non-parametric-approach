#working directory
setwd("~/Università/magistrale/Primo anno/Semestre 1.2/applied statistics/Progetto terremoti/prove diagnostics")

#import dataset
library(car)
library(plot3D)
library(rgl)
library(readr)
library(lme4)
library(nlme)
library(lattice)
library(influence.ME)
library(HLMdiag)
DATA<- read.csv("FAS_002_11_NoOutlier.csv")
head(DATA)
dim(DATA)
colnames(DATA)

# DATA<-DATA[sample(dim(DATA)[1],),]

#computing LogY
attach(DATA)
GeomMean=sqrt(FAS_EW*FAS_NS)
LogY=log10(GeomMean)     #already define the log for comfort
DATA<-cbind(DATA,GeomMean,LogY)
remove(GeomMean)
remove(LogY)

#to define the model, consider only the values for which Usable VS30 is defined
detach(DATA)
DATA<-DATA[-which(is.na(DATA$UsableVS30)),]
attach(DATA)


library(minpack.lm)
MODELLO<-nlsLM(LogY~const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
               + a2*ifelse((ML- mh)>0,ML-mh,0) 
               +c1*ML*log10(sqrt(Dipo^2+h^2))
               +kappa*log10(sqrt(Dipo^2+h^2))  ## kappa = -mref*c1 +c2
               +c3* sqrt(Dipo^2+h^2)
               +d1*log10(UsableVS30/800), 
               start=list(const=1.6,a1=1, a2=1.3, mh=4, h=3,c1=-0.1,c3=0.003,d1=-0.3,kappa=-4), control = nls.lm.control(maxiter = 1000),
)
MODELLO
summary(MODELLO)$coeff

c1<-summary(MODELLO)$coeff[6,1]
c3<-summary(MODELLO)$coeff[7,1]
h<-summary(MODELLO)$coeff[5,1]
mh<-summary(MODELLO)$coeff[4,1]
a1<-summary(MODELLO)$coeff[2,1]
a2<-summary(MODELLO)$coeff[3,1]
const<-summary(MODELLO)$coeff[1,1]
d1<-summary(MODELLO)$coeff[8,1]
kappa<-summary(MODELLO)$coeff[9,1]

kapparelation<-function(mref){kappa+c1*mref}
plot(min(ML):max(ML),kapparelation(min(ML):max(ML)))
lines(min(ML):max(ML),kapparelation(min(ML):max(ML)))
mref=4
c2<-kapparelation(mref)

coeff.estimates<-list(linear=c(const,a1,a2,c1,c2,c3,d1),nonlinear=c(mh,mref,h))

#clustering
detach(DATA)
de=dist(DATA$EvDpt, method = "euclidean")
eucl.avg=hclust(de, method = "average")
cluster.ea <- cutree(eucl.avg, k=2)  # euclidean-average
open3d()
Q<-cbind(EvLon=DATA$EvLon,EvLat=DATA$EvLat,EvDpt=DATA$EvDpt)
Q<-as.data.frame(Q)
plot3d(Q, size=3, col=cluster.ea+1, aspect = F)

#define groups
ST.fact<-as.factor(DATA$Stat)
Cluster.fact<-as.factor(cluster.ea-1)

#set nonlinear coefficients
mh     <-coeff.estimates[[2]][1]
mref   <-coeff.estimates[[2]][2]
h      <-coeff.estimates[[2]][3]

attach(DATA)
#define DESIGN MATRIX
Z<-data.frame(X0=rep(1,dim(DATA)[1]),X1=ifelse((ML- mh)<=0,ML-mh,0),X2=ifelse((ML- mh)>0,ML-mh,0),
              X3=(ML-mref)*log10(sqrt(Dipo^2+h^2)),X4=log10(sqrt(Dipo^2+h^2)),X5=sqrt(Dipo^2+h^2),
              X6=log10(UsableVS30/800))
attach(Z)
X<-cbind(Z,EV.fact,ST.fact,Cluster.fact,LogY)
LME.fit.cluster<-lmer(LogY~ X1+X2+X3+X4+X5+X6+Cluster.fact+(1|ST.fact),data=X)
summary(LME.fit.cluster)
brief(LME.fit.cluster)

#diagnostics
x11()
plot(LME.fit.cluster) 
x11()
qqnorm(residuals(LME.fit.cluster))
qqline(residuals(LME.fit.cluster))  #amazing
ResStat<-HLMresid(object=LME.fit.cluster,level="ST.fact")
res=residuals(LME.fit.cluster) # equal to res<-HLMresid(object=LME.fit.cluster,level=1)
lev<-hatvalues(LME.fit.cluster)
rstud=rstudent(LME.fit.cluster)     #residui studentizzati
# shapiro.test(res)
i1 <- which(Cluster.fact==0) #10047 surface earthquake
i2 <- which(Cluster.fact==1) #776 deep earthquake
