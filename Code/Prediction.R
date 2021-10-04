######## PREDICTION ########

setwd("C:/Users/aughi/Desktop/NONPARAMETRIC STAt/Project")
set.seed(210197)

## Build the classical model
library(readr)
library(lmerTest)

DATA<- read.csv("FAS_000_50_NoOutlier.csv")
load("C:/Users/aughi/Desktop/Progetto Stat/Dataset_Southern_Italy/Dataset_No_Outlier/coeff_estimates_000.50HZ.RData")

DATA<-DATA[-which(is.na(DATA$UsableVS30)),]
EV.fact<-as.factor(DATA$ID)
ST.fact<-as.factor(DATA$Stat)

mh     <-coeff.estimates[[2]][1]
mref   <-coeff.estimates[[2]][2]
h      <-coeff.estimates[[2]][3]

#MODEL:independent  mixed effects due to station and event
# LogY =  const+ a1*ifelse((ML- mh)<=0,ML-mh,0)
#        + a2*ifelse((ML- mh)>0,ML-mh,0) 
#        +(c1*(ML-mref)+c2)*log10(sqrt(Dipo^2+h^2))
#        +c3* sqrt(Dipo^2+h^2)
#        +d1*log10(UsableVS30/800)
#        +randomeffectdependingonST

#load packages
library(lme4)
library(nlme)
attach(DATA)

#define DESIGN MATRIX
Z<-data.frame(X0=rep(1,dim(DATA)[1]),X1=ifelse((ML- mh)<=0,ML-mh,0),X2=ifelse((ML- mh)>0,ML-mh,0),
              X3=(ML-mref)*log10(sqrt(Dipo^2+h^2)),X4=log10(sqrt(Dipo^2+h^2)),X5=sqrt(Dipo^2+h^2),
              X6=log10(UsableVS30/800))
attach(Z)

LME.fit<-lme4::lmer(DATA$LogY~ X1+X2+X3+X4+X5+X6+(1|ST.fact), 
              data=cbind(Z,EV.fact,ST.fact,DATA$LogY))
summary(LME.fit)

plot(LME.fit,xlab='Fitted',ylab='Residuals', main='Residuals Vs Fitted - Classic Model') #omoschedastici
qqnorm(residuals(LME.fit))
qqline(residuals(LME.fit),col='red') #Fat Right Tail

#plot on magnitude dipo grid
x.grid<-seq(range(DATA$ML)[1], range(DATA$ML)[2],length.out = 100)
y.grid<-seq(range(DATA$Dipo)[1], range(DATA$Dipo)[2],length.out = 100)
grid<-expand.grid(x.grid, y.grid)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
pred_vis<-matrix(data=NA,nrow=100, ncol=100)
for(i in 1:100){
  Z_vis<-data.frame(X0=x.grid[i],X1=ifelse((x.grid[i]- mh)<=0,x.grid[i]-mh,0),X2=ifelse((x.grid[i] - mh)>0,x.grid[i]-mh,0),
                   X3=(x.grid[i]-mref)*log10(sqrt(y.grid^2+h^2)),X4=log10(sqrt(y.grid^2+h^2)),X5=sqrt(y.grid^2+h^2),
                   X6=log10(mean(DATA$UsableVS30)/800),ST.fact=factor(getmode(DATA$Stat)))
  pred_vis[i,]<-predict(LME.fit,newdata=Z_vis, type='response')
}

library(rgl)
persp3d(x.grid,y.grid,pred_vis,col='dodgerblue2',xlab='Magnitude', ylab='Dipo', zlab='Regression',main='Magnitude-Hypocentral Distance Regression surface')
points3d(DATA$ML,DATA$Dipo,DATA$LogY, col='darkgrey', add=T)

library(plot3D)
persp3D(x.grid,y.grid,pred_vis,xlab='Magnitude', ylab='Dipo', zlab='Regression',main='Classic',theta=110,phi=10)
points3D(DATA$ML,DATA$Dipo,DATA$LogY, col='darkgrey',cex=0.02, add=T)

## BUild Nonparametric model
library(MASS)
library(rgl)
library(DepthProc)
library(gam)
library(plot3D)

DATA_NP<-read.csv("FAS_000_50_NoOutlier.csv", header=TRUE, sep=",")
Station=factor(DATA_NP$Stat)
attach(DATA_NP)

best_gam=mgcv::gam(LogY ~ s(Dipo,bs='cr')+s(ML,bs='cr')+s(EvDpt,bs='cr')+ s(Dipo,ML,bs='tp')+s(Station,bs="re")) 
summary(best_gam)

plot(best_gam$fitted.values,best_gam$residuals,xlab='Fitted', ylab='Residuals', col='dodgerblue', main='Residuals Vs Fitted - NP model')
abline(h=0.0, col='black')

qqnorm(best_gam$residuals)
qqline(best_gam$residuals,col='red') # same fat right tail

ks.test(x=best_gam$residuals,y='pnorm',alternative='two.sided') ###non-gaussian


##plot on magnitude dipo grid

pred_np_vis<-matrix(data=NA,nrow=100, ncol=100)
for(i in 1:100){
  Z_vis_np<-data.frame( ML=x.grid[i],Dipo=y.grid, EvDpt=mean(DATA$EvDpt), Station=factor(getmode(DATA$Stat)))
  pred_np_vis[i,]<-predict(best_gam,newdata=Z_vis_np)
}

library(rgl)
persp3d(x.grid,y.grid,pred_np_vis,sub='Nonparametric Model',col='green',xlab='Magnitude', ylab='Dipo', zlab='Regression',main='Magnitude-Hypocentral Distance Regression surface')
points3d(DATA$ML,DATA$Dipo,DATA$LogY, col='darkgrey', add=T)

library(plot3D)
par(mfrow=c(1,2))
persp3D(x.grid,y.grid,pred_vis,xlab='Magnitude', ylab='Dipo', zlab='Regression',main='Classic',theta=110,phi=10)
points3D(DATA$ML,DATA$Dipo,DATA$LogY, col='darkgrey',cex=0.02, add=T)
persp3D(x.grid,y.grid,pred_np_vis,sub='Nonparametric Model',xlab='Magnitude', ylab='Dipo', zlab='Regression',main='Non parametric',theta=110,phi=10)
points3D(DATA$ML,DATA$Dipo,DATA$LogY, col='darkgrey',cex=0.02, add=T)


## "Spatial" distribution of residuals comparison
par(mfrow=c(1,2))
scatter3D(DATA$ML, DATA$Dipo, DATA$LogY, pch = 18, cex = 0.01, colvar= residuals(LME.fit),
          theta = -40, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.grid, y = y.grid, z = pred_vis, col='darkgrey',
                      facets = NA, fit1 = predict(LME.fit)), main = "Classical")
scatter3D(DATA_NP$ML, DATA_NP$Dipo, DATA_NP$LogY, pch = 18, cex = 0.01, colvar = residuals(best_gam),
          theta = -40, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.grid, y = y.grid, z = pred_np_vis,col='darkgrey',  
                      facets = NA, fit1 = predict(best_gam)), main = "Nonparametric")
library(scatterplot3d)
scatterplot3d(DATA$ML, DATA$Dipo,residuals(LME.fit),pch='',type='h',color = as.numeric(cut(residuals(LME.fit),5)),
              xlab='ML', ylab='Dipo',zlab='Residuals', main='Classic model')


scatterplot3d(DATA_NP$ML, DATA_NP$Dipo,residuals(best_gam),pch='',type='h',color = as.numeric(cut(residuals(best_gam),5)),
              xlab='ML', ylab='Dipo',zlab='Residuals', main='Np model')
par(mfrow=c(1,2))
scatterplot3d(DATA$Dipo,DATA$ML,residuals(LME.fit),pch='',type='h',color = as.numeric(cut(residuals(LME.fit),5)),
              xlab='Dipo', ylab='ML',zlab='Residuals', main='Classic model')
scatterplot3d(DATA_NP$Dipo,DATA_NP$ML,residuals(best_gam),pch='',type='h',color = as.numeric(cut(residuals(best_gam),5)),
              xlab='Dipo', ylab='ML',zlab='Residuals', main='NP model')

## Two models comparison (fitting point of view)
#Considered Criteria: Aic, Bic,Aicc
print("AIC classical")
AIC(LME.fit)
print("AIC Nonparam")
AIC(best_gam)

print("BIC classical")
BIC(LME.fit)
print("BIC Nonparam")
BIC(best_gam)

library(qpcR)
print("AICc classical")
AICc(LME.fit)
print("AICc Nonparam")
AICc(best_gam)


## Prediction
#  Note that we cannot rely on gaussianity of residuals
library(lme4)

## load the test set:

test_set<-read.csv("FAS_000_test_NoOutlier.csv")
test_set<-test_set[-which(is.na(test_set$UsableVS30)),]
len<-dim(test_set)[1]
index<-sample(1:len,30)
test_set<-test_set[index,]
head(test_set)
actual<-test_set$LogY

attach(test_set)

#classic
Z_test<-data.frame(X0=rep(1,dim(test_set)[1]),X1=ifelse((ML- mh)<=0,ML-mh,0),X2=ifelse((ML- mh)>0,ML-mh,0),
              X3=(ML-mref)*log10(sqrt(Dipo^2+h^2)),X4=log10(sqrt(Dipo^2+h^2)),X5=sqrt(Dipo^2+h^2),
              X6=log10(UsableVS30/800),ST.fact=factor(test_set$Stat))

pred_class<-predict(LME.fit,newdata=Z_test, type='response')


#nonparam
Dipo_test=test_set$Dipo
ML_test  =test_set$ML
EvDpt_test=test_set$EvDpt
Station_test=factor(test_set$Stat)
newdata=data.frame(Dipo=Dipo_test, ML=ML_test, EvDpt=EvDpt_test, Station=Station_test)
pred_nonparam<-predict(best_gam,newdata=newdata,type='response')


#reverse percentile intervals
DATA_NP<-read.csv("FAS_000_50_NoOutlier.csv", header=TRUE, sep=",")
Station=factor(DATA_NP$Stat)
attach(DATA_NP)
library(pbapply)
set.seed(1234)
B<-1000
fitted.obs<-best_gam$fitted.values
res.obs <-best_gam$residuals

wrapper=function(){
  response.b=fitted.obs + sample(res.obs, replace=T)
  gam.b=mgcv::gam(response.b ~ s(Dipo,bs='cr')+s(ML,bs='cr')+s(EvDpt,bs='cr')+s(Dipo, ML, bs='tp')+s(Station,bs="re"))
  pred.b<-predict(gam.b, newdata=newdata, type='response')
}

pred.boot.L<-pbreplicate(B,wrapper(),simplify='vector')
load("C:/Users/aughi/Desktop/NONPARAMETRIC STAt/Project/re_perc_intervals_util.Rdata")
alpha<-0.05
r.q.L<- quantile(pred.boot.L, 1-alpha/2)
l.q.L<-quantile(pred.boot.L, alpha/2)
pred_nonparam<-as.vector(pred_nonparam)
CI<-cbind(pred_nonparam-(r.q.L-pred_nonparam), pred_nonparam, pred_nonparam-(l.q.L-pred_nonparam))

layout(1)
plot(CI[,2], col='green', cex=1.5, ylim=c(-8,1),pch=2, main='Classic Vs Nonparametric', xlab='', ylab='LogY')
points(CI[,1], col='blue', cex=1.7, pch='-')
points(CI[,3], col='blue', cex=1.7, pch='-')
segments(1:30,CI[,1],1:30,CI[,3], lwd=1)
points(pred_class, col='red', cex=1,pch=15)
points(actual, col='black', cex=1, pch=19)
x1<-1:30-0.2
x2<-1:30+0.2
y1<--7.99+abs(actual-pred_nonparam)
y2<--7.99+abs(actual-pred_class)
x.g<-rep(NA,60)
y.g<-rep(NA,60)
ind<-seq(1,60,by=2)
for  (i in 1:30){
  x.g[ind[i]]<-x1[i]
  x.g[ind[i]+1]<-x2[i]
  y.g[ind[i]]<-y1[i]
  y.g[ind[i]+1]<-y2[i]
}
library(plotfunctions)
add_bars(x.g,y.g,y0=-7.99,width=0.1,col=c('green','red'))
legend('topleft',c("nonparam", "classic", "true"), col=c('green','red','black'), pch=c(2,15,19))

RMSE <- function(error) { sqrt(mean(error^2)) }
RMSE_np<-round(RMSE(actual-pred_nonparam), digits = 4)
RMSE_cl<-round(RMSE(actual-pred_class), digits = 4)
par(mfrow=c(1,2))
hist(abs(actual-pred_nonparam),xlab = 'Error' ,col='darkgreen',main='Absolute errors - NP',breaks=30)
text(0.3,6,labels = paste("RMSE: ",RMSE_np))
hist(abs(actual-pred_class),xlab='Error', col='red',main='Absolute errors - Classic',breaks=30)
text(x=0.4,y=4.3,labels = paste("RMSE: ",RMSE_cl))

#GET FINAL RMSE ESTIMATE: 10-FOLD CV
set.seed(210197)
n<-dim(DATA)[1]
index<-sample(1:n, replace=F)
RMSE_NP_cv<-rep(NA,10)
RMSE_CL_cv<-rep(NA,10)
for(i in 1:10){
  test_<-((i-1)*999+1):((i-1)*999+999)
  test_ind<-index[test_]
  train_ind<-index[-test_]
  testData<-DATA[test_ind,]
  trainData<-DATA[-test_ind,]
  
  ## Model training
  Z<-data.frame(X0=rep(1,dim(trainData)[1]),X1=ifelse((trainData$ML- mh)<=0,trainData$ML-mh,0),X2=ifelse((trainData$ML- mh)>0,trainData$ML-mh,0),
                          X3=(trainData$ML-mref)*log10(sqrt(trainData$Dipo^2+h^2)),X4=log10(sqrt(trainData$Dipo^2+h^2)),X5=sqrt(trainData$Dipo^2+h^2),
                          X6=log10(trainData$UsableVS30/800))
  ST.fact<-as.factor(trainData$Stat)
  
  OldModel<-lme4::lmer(trainData$LogY~ X1+X2+X3+X4+X5+X6+(1|ST.fact), 
                      data=cbind(Z,as.factor(trainData$ID),ST.fact,trainData$LogY))
  LogY=trainData$LogY
  Dipo=trainData$Dipo
  ML  =trainData$ML
  EvDpt=trainData$EvDpt
  Station<-ST.fact
  NpModel <-mgcv::gam(LogY ~ s(Dipo,bs='cr')+s(ML,bs='cr')+s(EvDpt,bs='cr')+s(Dipo,ML, bs='tp')+s(Station,bs="re"))
  
  ## Predictions
  ac<-testData$LogY
  Z_t<-data.frame(X0=rep(1,dim(testData)[1]),X1=ifelse((testData$ML- mh)<=0,testData$ML-mh,0),X2=ifelse((testData$ML- mh)>0,testData$ML-mh,0),
                     X3=(testData$ML-mref)*log10(sqrt(testData$Dipo^2+h^2)),X4=log10(sqrt(testData$Dipo^2+h^2)),X5=sqrt(testData$Dipo^2+h^2),
                     X6=log10(testData$UsableVS30/800),ST.fact=factor(testData$Stat))
  
  pred_class<-predict(OldModel,newdata=Z_t, type='response')
  
  newdata1=data.frame(Dipo=testData$Dipo, ML=testData$ML, EvDpt=testData$EvDpt, Station=factor(testData$Stat))
  pred_nonparam<-predict(NpModel,newdata1,type='response')
  
  ##RMSE ESTIMATEhe
  RMSE_CL_cv[i]<-RMSE(ac-pred_class)
  RMSE_NP_cv[i]<-RMSE(ac-pred_nonparam)
}

RMSE_CL<-sum(RMSE_CL_cv)/10
RMSE_NP<-sum(RMSE_NP_cv)/10
print("CV estimates are: ")
print(paste("Classic Approach: ",RMSE_CL))
print(paste("Nonparam Approach: ",RMSE_NP))

