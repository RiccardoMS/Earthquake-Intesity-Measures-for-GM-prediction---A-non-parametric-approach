#working directory
setwd("C:/Users/aughi/Desktop/Progetto Stat/Dataset_Southern_Italy/Dataset_No_Outlier")

#import dataset
library(readr)
DATA<- read.csv("FAS_000_50_NoOutlier.csv")
head(DATA)
dim(DATA)
colnames(DATA)
attach(DATA)

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

coeff.estimates.000_50HZ<-list(linear=c(const,a1,a2,c1,c2,c3,d1),nonlinear=c(mh,mref,h))
save(coeff.estimates, file="coeff_estimates_000.50HZ.RData")

##note that coefficients are saved in a LIST
## TO ACCESS:
#  load the object(double click on his name in the working directory
#                 or load("percorso/file/coeff.estimates.000_50HZ.RData") in the script
#  coeff.estimates.000_50HZ[[1]][1 to 7] for the linear coefficients    (const-a1-a2-c1-c2-c3-d1)
#  coeff.estimates.000_50HZ[[2]][1 to 3] for the non linear coefficients (mh--mref--h)