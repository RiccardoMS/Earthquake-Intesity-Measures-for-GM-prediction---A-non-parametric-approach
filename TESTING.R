##############################      SCRIPT 1      ######################################
#working directory
setwd("C:/Users/aughi/Desktop/NONPARAMETRIC STAt/Project")

#import dataset
DATA<-read.csv("FAS_000_50.csv", header=TRUE, sep=",")
head(DATA)
dim(DATA)
colnames(DATA)


#define GeomMean of FAS_EW, FAS_NS
attach(DATA)
GeomMean=sqrt(DATA$FAS_EW*DATA$FAS_NS)
LogY=log10(GeomMean)     
DATA<-cbind(DATA,GeomMean,LogY)
detach(DATA)
remove(GeomMean)
remove(LogY)
attach(DATA)

GeomMean.mean=mean(GeomMean)                  #variable of interest but problem in distribution
GeomMean.sd=sd(GeomMean)
par(mfrow=c(1,2))
plot(GeomMean,main='Log of Geom Mean')  
abline(a=GeomMean.mean,b=0,col='green')   
abline(a=GeomMean.mean+GeomMean.sd,b=0,col='red')
hist(GeomMean,main = 'Distribution of GeomMean',prob=T)
qqnorm(GeomMean)                          
qqline(GeomMean)
boxplot(GeomMean) 

LogY.mean=mean(LogY)                          #better distribution of data, but assumption of
LogY.sd=sd(LogY)                              #normality of the response is not completely true
par(mfrow=c(1,2))
plot(LogY,main='Log of Geom Mean')  
abline(a=LogY.mean,b=0,col='green')   
abline(a=LogY.mean+LogY.sd,b=0,col='red')
hist(LogY,main = 'Distribution of Log of GeomMean',prob=T)
qqnorm(LogY)                          
qqline(LogY)
boxplot(LogY)  

B<-1000
shapiro.pvalue<-numeric(B)
for (b in 1:B){
  Yobs<-sample(LogY,5000)
  shapiro.pvalue[b]<-shapiro.test(Yobs)$p.value
}
sum(shapiro.pvalue>=0.1)/B                     #refuse normality



#DIPO
Dipo.mean=mean(Dipo)
Dipo.sd=sd(Dipo)
par(mfrow=c(1,2))
plot(Dipo,main='occurencies of DIPO')  
abline(a=Dipo.mean,b=0,col='green')   
abline(a=Dipo.mean+Dipo.sd,b=0,col='red')
hist(Dipo,main = 'distribution of occurencies',prob=T)
qqnorm(Dipo)                          
qqline(Dipo)
boxplot(Dipo) 

table(Dipo>200)                            

detach(DATA)
DATALONGDISTANCE<-DATA[which(DATA$Dipo>200),]
DATA<-DATA[-which(DATA$Dipo>200),]
attach(DATA)

#define Usable VS_30            
UsableVS30=ifelse(is.na(DATA$VS30),DATA$VS30_WA,DATA$VS30)
DATA<-cbind(DATA,UsableVS30)
remove(UsableVS30)

table(is.na(DATA$UsableVS30))

DATANOVS30<-DATA[which(is.na(DATA$UsableVS30)),]
DATA<-DATA[-which(is.na(DATA$UsableVS30)),]
DATA<-DATA[-c(2701,2366),]
head(DATA)
head(DATALONGDISTANCE)
head(DATANOVS30)

DATAID<-DATA[,1]
DATASTAT<-DATA[,7]
DATAEC8<-DATA[,17]
DATA<-DATA[,-c(1,3,4,6,7,9,10,11,13,14,15,16,17,18,19,20)]
DATALONGDISTANCE<-DATALONGDISTANCE[,-c(1,3,4,6,7,9,10,11,13,14,15,16,17,18,19,20)]
DATANOVS30<-DATANOVS30[,-c(1,3,4,6,7,9,10,11,13,14,15,16,17,18,19,20)]


############################# RANKINGS & DEPTHS ##########################################
library(depth)
library(DepthProc)

prova=depth(DATA,method='Tukey')
hist(prova)

depthMedian(DATA,depth_params = list(method='Tukey')) # ALL VARIABLES INTO ACCOUNT

depthContour(DATA[,c(1,5)],depth_params = list(method='Tukey')) #workonly on bivariates, i dont think it makes any sense to do so

bgplot=aplpack::bagplot.pairs(DATA)

bgplotML_Mstat=aplpack::bagplot(DATA$ML,DATA$MStat)
bgplotML_Mstat$pxy.outlier

bgplotML_Logy=aplpack::bagplot(DATA$ML,DATA$LogY)
bgplotML_Logy$pxy.outlier

bgplotEvDpt_Logy=aplpack::bagplot(DATA$EvDpt,DATA$LogY,xlab='Event Depth', ylab='LogY')
bgplotEvDpt_Logy$pxy.outlier

bgplotMStat_Logy=aplpack::bagplot(DATA$MStat,DATA$LogY)
bgplotMStat_Logy$pxy.outlier

bgplotDipo_Logy=aplpack::bagplot(DATA$Dipo,DATA$LogY, xlab='Hypocentral Distance', ylab='LogY')
bgplotDipo_Logy$pxy.outlier

bgplotUsableVS30_Logy=aplpack::bagplot(DATA$UsableVS30,DATA$LogY)
bgplotUsableVS30_Logy$pxy.outlier


#################################### TESTING ###########################################
set.seed(210197)
detach(DATA)
attach(DATA)
attach(DATAEC8)
## test on LogY distribution in LongDistance Data

depthMedian(DATALONGDISTANCE,depth_params = list(method='Tukey'))

# H0: DATA =d= LONGDISTANCEDATA vs H1: DATA=!d!= LONGDISTANCEDATA

T0<-abs(mean(DATA$LogY)-mean(DATALONGDISTANCE$LogY))
DataPooled<-c(DATA$LogY,DATALONGDISTANCE$LogY)
n<-length(DataPooled)
B<-10000
tstat<-numeric(B)
for(perm in 1:B){
  permutation<-sample(1:n)
  x_perm<-DataPooled[permutation]
  x1_perm<-x_perm[1:length(DATA$LogY)]
  x2_perm<-x_perm[(length(DATA$LogY)+1):n]
  tstat[perm]<-abs(mean(x1_perm)-mean(x2_perm))
}
hist(tstat,xlim=range(c(tstat,T0)),breaks = 20)
abline(v=T0,col='red')
plot(ecdf(tstat))
pvalue<-sum(tstat>=T0)/B
pvalue

#what if we include Dipo?
t1<-DATA[,-c(2,6)]
t2<-DATALONGDISTANCE[,-2]
t1.mean<-colMeans(t1)
t2.mean<-colMeans(t2)
n1<-dim(t1)[1]
n2<-dim(t2)[2]
T0<-as.numeric((t1.mean-t2.mean)%*%(t1.mean-t2.mean))
DataPooled<-rbind(t1,t2)
n<-n1+n2
B<-10000
tstat<-numeric(B)
for(perm in 1:B){
  permutation<-sample(n)
  x_perm<-DataPooled[permutation,]
  x1_perm<-x_perm[1:n1,]
  x2_perm<-x_perm[(n1+1):n,]
  x1.mean<-colMeans(x1_perm)
  x2.mean<-colMeans(x2_perm)
  tstat[perm]<-as.numeric((x1.mean-x2.mean)%*%(x1.mean-x2.mean))
}
hist(tstat,xlim=range(c(tstat,T0)),breaks = 20)
abline(v=T0,col='red')
plot(ecdf(tstat))
pvalue<-sum(tstat>=T0)/B
pvalue
## dramatic change: actually make sense to think of two different distribution
##                  if the aim is to train a model







#Hypothesis: MStat=f(ML)--> if so, right to not consider one of them
result<-lm(DATA$MStat ~ DATA$ML)
summary(result)
qqnorm(result$residuals)

T01<-summary(result)$f[1]
B<-10000
TH01<-numeric(B)
for (i in 1:B){
  perm<-sample(length(DATA$MStat))
  yperm<-DATA$MStat[perm]
  result_perm<-lm(yperm~DATA$ML)
  TH01[i]<-summary(result_perm)$f[1]
}
pvalue<-sum(TH01>=T01)/B
hist(TH01)
abline(v=T01)             
plot(ecdf(TH01))

## verify effect of ec8
layout(1)
boxplot(LogY~DATAEC8, col=rainbow(8),main='Effect of EC8 on LogY', xlab='EC8')
#hypothesys of normality in groups dissatisfied: permutational anova
g<-nlevels(DATAEC8)
n<-length(DATAEC8)
fit <- aov(LogY ~ DATAEC8)
summary(fit)

T0 <- summary(fit)[[1]][1,4]
T0

B <- 1000 # Number of permutations
T_stat <- numeric(B) 

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  LogY_perm<- LogY[permutation]
  fit_perm <- aov(LogY_perm ~ DATAEC8)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}


hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-1,20))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val

#test A against B

ind1<-which(DATAEC8=='A')
ind2<-which(DATAEC8=='B')
data1<-cbind(LogY[ind1], DATAEC8[ind1])
data2<-cbind(LogY[ind2], DATAEC8[ind2])
data<-rbind(data1,data2)
data<-as.data.frame(data)
n<-length(data[,1])
fit <- aov(data[,1] ~ data[,2])
summary(fit)

T0 <- summary(fit)[[1]][1,4]
T0

B <- 1000 # Number of permutations
T_stat <- numeric(B) 

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  LogY_perm<- data[permutation,1]
  fit_perm <- aov(LogY_perm ~ data[,2])
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}


hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-1,20))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val

## Astar against Bstar

ind1<-which(DATAEC8=='A*')
ind2<-which(DATAEC8=='B*')
data1<-cbind(LogY[ind1], DATAEC8[ind1])
data2<-cbind(LogY[ind2], DATAEC8[ind2])
data<-rbind(data1,data2)
data<-as.data.frame(data)
n<-length(data[,1])
fit <- aov(data[,1] ~ data[,2])
summary(fit)

T0 <- summary(fit)[[1]][1,4]
T0

B <- 1000 # Number of permutations
T_stat <- numeric(B) 

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  LogY_perm<- data[permutation,1]
  fit_perm <- aov(LogY_perm ~ data[,2])
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}


hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-1,20))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val

#test A against E

ind1<-which(DATAEC8=='A')
ind2<-which(DATAEC8=='E')
data1<-cbind(LogY[ind1], DATAEC8[ind1])
data2<-cbind(LogY[ind2], DATAEC8[ind2])
data<-rbind(data1,data2)
data<-as.data.frame(data)
n<-length(data[,1])
fit <- aov(data[,1] ~ data[,2])
summary(fit)

T0 <- summary(fit)[[1]][1,4]
T0

B <- 1000 # Number of permutations
T_stat <- numeric(B) 

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  LogY_perm<- data[permutation,1]
  fit_perm <- aov(LogY_perm ~ data[,2])
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}


hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-1,20))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val


##test diffusive hypothesis for UsableVS30
bgplotUsableVS30_Logy=aplpack::bagplot(DATA$UsableVS30,DATA$LogY)


result<-lm(LogY~UsableVS30)

qqnorm(result$residuals)
qqline(result$residuals)

n<-length(LogY)
T0<-summary(result)$f[1]
T0

B<-1000
T_stat<-numeric(B)

for (perm in 1:B){
  permutazione <- sample(n)
  LogY_perm<-LogY[permutazione]
  T_stat[perm] <- summary(lm(LogY_perm~UsableVS30))$f[1]
}

hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat),xlim=c(-1,20))
abline(v=T0,col=3,lwd=4)

p_val <- sum(T_stat>=T0)/B
p_val


### evdpt
bgplotUsableVS30_Logy=aplpack::bagplot(DATA$EvDpt,DATA$LogY)


result<-lm(DATA$LogY ~ DATA$EvDpt)
summary(result)
qqnorm(result$residuals)

T01<-summary(result)$f[1]
B<-10000
TH01<-numeric(B)
for (i in 1:B){
  perm<-sample(length(DATA$LogY))
  yperm<-DATA$LogY[perm]
  result_perm<-lm(yperm~DATA$EvDpt)
  TH01[i]<-summary(result_perm)$f[1]
}
pvalue<-sum(TH01>=T01)/B
hist(TH01)
abline(v=T01)                 


plot(ecdf(TH01))


##multiple regression

result <- lm(LogY ~ MStat + Dipo + UsableVS30+ EvDpt)
summary(result)
qqnorm(result$residuals)


# H0: beta3 = 0
# test statistic
summary(result)$coefficients
T0_x3 <- abs(summary(result)$coefficients[4,3])
T0_x3
# permutations
# residuals of the reduced model
# reduced model:
# Y = beta0 + beta1*x1 + beta2*x2
regr.H03 <- lm(LogY ~ MStat + Dipo+ EvDpt)
residui.H03 <- regr.H03$residuals
residui.H03.perm <- residui.H03[permutazione]
# permuted y:
Y.perm.H03 <- regr.H03$fitted + residui.H03.perm


# H0: beta4 = 0
# test statistic
summary(result)$coefficients
T0_x4 <- abs(summary(result)$coefficients[5,3])
T0_x4
# permutations
# residuals of the reduced model
# reduced model:
regr.H04 <- lm(LogY ~ MStat + Dipo+ UsableVS30)
residui.H04 <- regr.H03$residuals
residui.H04.perm <- residui.H04[permutazione]
# permuted y:
Y.perm.H04 <- regr.H04$fitted + residui.H04.perm


#pvalues

B <- 1000
T_H03 <- T_H04 <- numeric(B)

for(perm in 1:B){
  permutazione <- sample(n)
  
  residui.H03.perm <- residui.H03[permutazione]
  Y.perm.H03 <- regr.H03$fitted + residui.H03.perm
  T_H03[perm] <- abs(summary(lm(Y.perm.H03 ~ MStat + Dipo + UsableVS30 +EvDpt))$coefficients[4,3])
  
  residui.H04.perm <- residui.H04[permutazione]
  Y.perm.H04 <- regr.H04$fitted + residui.H04.perm
  T_H04[perm] <- abs(summary(lm(Y.perm.H04 ~ MStat + Dipo + UsableVS30 +EvDpt))$coefficients[5,3])
  
}

sum(T_H03>=T0_x3)/B
sum(T_H04>=T0_x4)/B