#working directory
setwd("C:/Users/aughi/Desktop/NONPARAMETRIC STAt/Project")

#import dataset
DATA<-read.csv("FAS_000_50.csv", header=TRUE, sep=",")
head(DATA)
dim(DATA)
colnames(DATA)

##note that:
#ID         :=        identifies events,useful later for grouping
#MStat      :=        local magnitude
#EvDPt      :=        depth of event in km [useful later for clustering?]
#EvLAt-EvLon:=        coordinates of event,useful only for mapping thedataset
#StLat-StLon:=        coordinates of station where the record was taken,see above
#Stat       :=        code identifying stations, useful later for grouping
#Net        :=        identifies network of the station, till now no practical usefulness known
#Depi       :=        epicentral distance in km (the0rically applying pithagora to Depi + EvDpt we obtain DIpo?)
#Dipo       :=        hipocentral distance in km
#Azi        :=        source-site azimut
#FAS_Z/NS/EW:=        component of FAS; vertical component Z is included but not interesting for the future model
#EC8        :=        cathegorical variable, identify the type of earthquake
#VS30       :=        parameter relative to crustal properties, high number of NA values
#VS30_WA    :=        proxy for VS30, SOME INSTANCES DOES NOT HAVE NEITHER VS30 NOR VS30_WA!

attach(DATA)

########################## MAPPA DATASET ##########################################
library(rworldmap) #requires to download the package



newmap <- getMap(resolution = "high")

plot(newmap,
     xlim = c(12.3, 18.3),
     ylim = c(36.8,42.2),
     asp = 1,
     main = 'Dataset Visualization'
)

points(EvLon,EvLat,type = 'p', col='red',cex=0.6 )
points(StLon,StLat,col='blue',pch=4,cex=0.6)
leg.txt=c("Events","Stations")
legend("topright",leg.txt,col = c("red","blue"),pch='o',title='Legend')

##############################  DATASET ANALISYS ######################################
#ID-STAT
EventGroup=factor(ID)   #593 events
StatGroup=factor(Stat)  #316 stations

#FAS
FAS_Z.mean=mean(FAS_Z)
FAS_Z.sd=sd(FAS_Z)
par(mfrow=c(1,2))
plot(FAS_Z,main='occurencies in DATA') #presence of a clear outlier
identify(FAS_Z)                        #point 2701
abline(a=FAS_Z.mean,b=0,col='green')
abline(a=FAS_Z.mean+FAS_Z.sd,b=0,col='red')
hist(FAS_Z,main = 'distribution of occurencies',prob=T)#very unbalanced distribution, try log10() transform
hist(log10(FAS_Z),main = 'distribution of occurencies',prob=T) #clearly better
qqnorm(log10(FAS_Z),ylim=c(-5,2))                              #gaussianity is not assured in tails
qqline(log10(FAS_Z))

FAS_EW.mean=mean(FAS_EW)
FAS_EW.sd=sd(FAS_EW)
par(mfrow=c(1,2))
plot(FAS_EW,main='occurencies in DATA') #presence of outlier
identify(FAS_EW)                        #same point as before
abline(a=FAS_EW.mean,b=0,col='green')
abline(a=FAS_EW.mean+FAS_EW.sd,b=0,col='red')
hist(FAS_EW,main = 'distribution of occurencies',prob=T) #unbalanced,see above
hist(log10(FAS_EW),main = 'distribution of occurencies',prob=T)#clearly better 
qqnorm(FAS_EW,ylim=c(-5,2))                     #also here, gaussianity not assured, very heavy tails
qqline(FAS_EW)

FAS_NS.mean=mean(FAS_NS)
FAS_NS.sd=sd(FAS_NS)
par(mfrow=c(1,2))
plot(FAS_NS,main='occurencies in DATA') #presence of outlier
identify(FAS_NS)                        #same point as before
abline(a=FAS_NS.mean,b=0,col='green')
abline(a=FAS_NS.mean+FAS_NS.sd,b=0,col='red')
hist(FAS_NS,main = 'distribution of occurencies',prob=T) #unbalanced,see above
hist(log10(FAS_NS),main = 'distribution of occurencies',prob=T)#clearly better 
qqnorm(log10(FAS_NS),ylim=c(-5,2))                     #also here, gaussianity not assured, very heavy tails
qqline(log10(FAS_NS))

#control on point DATA[2701]
DATA[2701,] #ML=5.07, DEPI==DEPO=5,3; big magnitude small distance?
range(ML)   #not so high with respect to others
range(Dipo) #not so small wrt to others
par(mfrow=c(1,2))
plot(ML,FAS_EW)
identify(ML,FAS_EW)   #the combining of the two effects can result in a abnormous value 
plot(Dipo,FAS_EW)     #we decide for now to remove it
identify(Dipo,FAS_EW)

detach(DATA)
DATA<-DATA[-2701,]
attach(DATA)

par(mfrow=c(2,3))
plot(FAS_Z,main='occurencies in DATA')
plot(FAS_EW,main='occurencies in DATA')
plot(FAS_NS,main='occurencies in DATA')
plot(log10(FAS_Z),main='occurencies in DATA')     #we actually decided to work on a log10 scale
plot(log10(FAS_EW),main='occurencies in DATA')    #for FAS variables
plot(log10(FAS_NS),main='occurencies in DATA')    #note it is actually sensed according to the fact 
                                                  #that the variable of interest in the model
                                                  #is defined as the log10(geometricmean)
                                                  


#ML
ML.mean=mean(ML)
ML.sd=sd(ML)
par(mfrow=c(1,2))
plot(ML,main='occurencies in DATA') #plot appears with less occurencies simply 'cause the values are 
                                    #precise till 2nd decimal,so they're repeated
abline(a=ML.mean,b=0,col='green')   #distribution appears good enough
abline(a=ML.mean+ML.sd,b=0,col='red')
hist(ML,main = 'distribution of occurencies',prob=T)
qqnorm(ML)                          #non normality: trying some transformation
qqline(ML)
boxplot(ML)

plot(log10(ML),main='occurencies in DATA') 
hist(log10(ML),main = 'distribution of occurencies',prob=T)
qqnorm(log10(ML))                          #actually the distribution is slightly better, not so much
qqline(log10(ML))                          #could think about using log10 o raw data
boxplot(log10(ML))                        #better behavoiur of head tail


#MStat
MStat.mean=mean(MStat)
MStat.sd=sd(MStat)
par(mfrow=c(1,2))
plot(MStat,main='occurencies in DATA')  
abline(a=MStat.mean,b=0,col='green')   
abline(a=MStat.mean+MStat.sd,b=0,col='red')
hist(MStat,main = 'distribution of occurencies',prob=T)
qqnorm(MStat)                          
qqline(MStat)
boxplot(MStat)

plot(log10(MStat),main='occurencies in DATA') 
hist(log10(MStat),main = 'distribution of occurencies',prob=T)
qqnorm(log10(MStat))                         
qqline(log10(MStat))                          
boxplot(log10(MStat))      #overall behaviour of MStat appears to be like ML, some differences
                           #ML should be the true (estimated) magnityde of an event, MStat how it was
                           #percepted at the station. How to use MStat?

#DIPO
Dipo.mean=mean(Dipo)
Dipo.sd=sd(Dipo)
par(mfrow=c(1,2))
plot(Dipo,main='occurencies in DATA')  
abline(a=Dipo.mean,b=0,col='green')   
abline(a=Dipo.mean+Dipo.sd,b=0,col='red')
hist(Dipo,main = 'distribution of occurencies',prob=T)
qqnorm(Dipo)                          
qqline(Dipo)
boxplot(Dipo)            #clear evidence of small number of occurencies with anomalous distribution


#DEPI
Depi.mean=mean(Depi)
Depi.sd=sd(Depi)
par(mfrow=c(1,2))
plot(Depi,main='occurencies in DATA')  
abline(a=Depi.mean,b=0,col='green')   
abline(a=Depi.mean+Depi.sd,b=0,col='red')
hist(Depi,main = 'distribution of occurencies',prob=T)
qqnorm(Depi)                          
qqline(Depi)
boxplot(Depi)            #behaviour similar to DIPO


par(mfrow=c(1,2))
plot(Dipo)                             #very small number of occurencies between 300-500?
plot(Depi)                             #same point as below
table(Dipo>200)
table(Depi>200)                        #ARE THE SAME? SEEMS SO
table((Depi>200)==(Dipo>200))          #It is
pie(table(Dipo>200),labels= list("under 200","over 200"),col=c("red","blue"),main="DIPO")
pie(table(Depi>200),labels= list("under 200","over 200"),col=c("red","blue"),main="DEPI")

#Due to the fact that the cardinality of the group 'HIGH DISTANCE'is only 16 over >11000,we remove them 
detach(DATA)
DATA<-DATA[-which(DATA$Dipo>200),]
attach(DATA)
#Azi
Azi.mean=mean(Azi)       
Azi.sd=sd(Azi)
par(mfrow=c(1,2))
plot(Azi,main='occurencies in DATA')  
abline(a=Azi.mean,b=0,col='green')   
abline(a=Azi.mean+Azi.sd,b=0,col='red')
hist(Azi,main = 'distribution of occurencies',prob=T)
qqnorm(Azi)                          
qqline(Azi)
boxplot(Azi)                   #in my opinion, completely unuseful



#EC8
barplot(table(EC8)/length(EC8),ylim=c(0,1),xaxt="n")
axis(1,1:8,levels(EC8))
pie(table(EC8))   #practically Astar and Bstar cover the 0.76 of the data, adding the NA the 85 pc
                 #explore the relationship with response, maybe the belonging to a special class
                  #|E|=10; |A|=272; |B|=728; |C|=352; has an influence on the response

#VS30/VS30_WA
table(is.na(VS30))  
levels(factor(VS30))  #only 1191 values available relative to 61 sites

table(is.na(VS30_WA)) #9991 values availables, 1160 are not
#relationship between the two? WA seems NOT to be a proxy for VS30:
campione<-DATA[-which(is.na(VS30)),]
plot(abs(campione$VS30-campione$VS30_WA))


#DEFINE THE RESPONSE: GEOMETRIC MEAN OF FAS_EW & FAS__NS
detach(DATA)
GeomMean=sqrt(DATA$FAS_EW*DATA$FAS_NS)
LogY=log10(GeomMean)     #already define the log for cpmfort
DATA<-cbind(DATA,GeomMean,LogY)
remove(GeomMean)
remove(LogY)
attach(DATA)

LogY.mean=mean(LogY)       
LogY.sd=sd(LogY)
par(mfrow=c(1,2))
plot(LogY,main='occurencies in DATA')  
abline(a=LogY.mean,b=0,col='green')   
abline(a=LogY.mean+LogY.sd,b=0,col='red')
hist(LogY,main = 'distribution of occurencies',prob=T)
qqnorm(LogY)                          
qqline(LogY)
boxplot(LogY)  

#TRANSFORM CATHEGORICAL VAR EC8 IN SEQUENCE OF DUMMY VARIABLES
EC8vec=as.vector(EC8)
DummyA=EC8vec
DummyAStar=EC8vec
DummyB=EC8vec
DummyBStar=EC8vec
DummyC=EC8vec
DummyCStar=EC8vec
DummyE=EC8vec
lungh=length(EC8vec)
for (i in (1:lungh)){
  if(EC8vec[i]=="A"){
    DummyA[i]=1
    }
  else
    DummyA[i]=0
  
  if(EC8vec[i]=="A*"){
    DummyAStar[i]=1
    }
  else
    DummyAStar[i]=0
  
  if(EC8vec[i]=="B"){
    DummyB[i]=1
    }
  else
    DummyB[i]=0
  
  if(EC8vec[i]=="B*"){
    DummyBStar[i]=1
    }
  else
    DummyBStar[i]=0
  
  if(EC8vec[i]=="C"){
    DummyC[i]=1
    }
  else
    DummyC[i]=0
  
  if(EC8vec[i]=="C*"){
    DummyCStar[i]=1
    }
  else
    DummyCStar[i]=0
  
  if(EC8vec[i]=="E"){
    DummyE[i]=1
    }
  else
    DummyE[i]=0
}

DATA<-cbind(DATA,DummyA,DummyAStar,DummyB,DummyBStar,DummyC,DummyCStar,DummyE)
remove(DummyA)
remove(DummyAStar)
remove(DummyB)
remove(DummyBStar)
remove(DummyC)
remove(DummyCStar)
remove(DummyE)

##DEFINE USABLE_VS30             (Assign 1500 to NaN occurencies, proviisonal)
detach(DATA)
UsableVS30=ifelse(is.na(DATA$VS30),DATA$VS30_WA,DATA$VS30)
DATA<-cbind(DATA,UsableVS30)
remove(UsableVS30)
attach(DATA)
table(is.na(UsableVS30))


###################################### MAHALANOBIS DISTANCE##########################################

#define quantitative DATASET (WE DON'T INCLUDE SOME VARIABLES: AZI is completely irrelevant,moreover ia funtion 
#Depi (or Dipo) and EvDpt, VS30 and VS30 WA admits a too large numbers of NA occurencies; all lats and longs)
library(eqs2lavaan)
quant.DATA<-DATA[,c(2,5,8,11,12,13,14,15,16,21)]
m.quant.DATA=colMeans(quant.DATA)
COV.quant.DATA=cov(quant.DATA)
plotCov(COV.quant.DATA)     #actually suggests to avoid considering AZi, both depi and Dipo
remove(quant.DATA)

#try to consider instead the dataset of interest for the future MEL model:
#ML,Dipo,Mstat,LogY,EvDpt

model.DATA<-DATA[,c(2,5,8,12,21)]
m.model.DATA=colMeans(model.DATA)
COV.model.DATA=cov(model.DATA)
plotCov(COV.model.DATA)
maha.dist.model=mahalanobis(model.DATA,m.model.DATA,COV.model.DATA)
plot(maha.dist.model)
identify(maha.dist.model)      

####################################### RELATIONSHIP BETWEEM VARIABLES ################################
#overview
pairs(model.DATA)     #actually ML appears to be directly correlated with MStat

#ML-MStat
res=ML-MStat
plot(res)
abline(a=mean(res),b=0,col='green')
abline(a=mean(res)+2.5*sd(res),b=0,col='red')
abline(a=mean(res)-2.5*sd(res),b=0,col='red')
problematic.points<-DATA[which((abs(res)>2.5*sd(res))==T),]
prob.res=problematic.points$ML-problematic.points$MStat
points(which((abs(res)>2.5*sd(res))==T),prob.res,col='red')    #USEFUL OR NOT?

#ML-logY (Dipo bins)
range(Dipo)
Dipo.lab=c("0-50","50-100","100-150")
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
#nota: una dipendenza lineare sembra giustificata. Due problemi:
#leggero gomito tra ml=5 e ml=5.5
#dispersione elevata per magnitudo basse nel gruppo a distanza ridotta 0-50km

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
plot(log10(Dipo),LogY, main='Attenuation of Responde with log10Dipo',xlab='Log10Dipo[km]',ylab='Log10 of GeomMean',col=ML.col)
legend('bottomright',ML.labels,col=col.span,pch='o')
#to evaluate:possible transofrmation of DIpo in log(DIpo)


#LogY-EC8
pie(table(EC8=="")) # 1070 instnaces are not classified
boxplot(LogY~EC8)   #belonging to class c cstar e suggests a possible cluster,
                    #not so many infos provided

#LogY-Vs30 (????)
plot(VS30,LogY)
plot(VS30_WA,LogY)  #no idea about how to interpret


#LogY=function(ML,Dipo)  (experimental)
library('plot3D')

#distribution of occurencies in relation to a ML-Dipo Grid
hist3D_fancy(ML,Dipo,colvar=LogY,breaks=20,phi=10, theta=30)
hist3D_fancy(ML,Dipo,colvar=LogY,breaks=20,phi = 10,theta=-30)
#prevalence of instances in which ml and dipo are attained, conrespondent to values of logy small (-3,-4)

scatter3D(ML,Dipo,LogY,bty="g",ticktype="detailed",theta=20,phi=15)
scatter3D(ML,Dipo,LogY,bty="g",ticktype="detailed",theta=40,phi=15)
scatter3D(ML,Dipo,LogY,bty="g",ticktype="detailed",theta=60,phi=15)
scatter3D(ML,Dipo,LogY,bty="g",ticktype="detailed",theta =120,phi=15)
scatter3D(ML,Dipo,LogY,bty="g",ticktype="detailed",theta=140,phi=15)
scatter3D(ML,Dipo,LogY,bty="g",ticktype="detailed",theta=160,phi=15)
#possible plane of regression? actually seems to be a curve at extreme points
#regression plane
# x, y, z variables
x <- ML
y <- Dipo
z <- LogY
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 0.25, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA), main = "RegPlane")
scatter3D(x, y, z, pch = 18, cex = 0.01, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "RegPlane")
scatter3D(x, y, z, pch = 18, cex = 0.01, 
          theta = 0, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "RegPlane")
summary(fit)
#DIpo as an ininfluent coefficient--> try some transformation
#model seems to not fit well particularly the instnaces with high values of LogY
plot(LogY,fit$residuals) #actually it is so

#try a ITA18-like model
x <- ML
y <- Dipo
z <- LogY
# Compute the linear regression (z = ax +cx*log10(y) + dy+ e)
fit1 <- lm(z ~ x + x*log10(y) + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit1, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fit1ted points for droplines to surface
fit1points <- predict(fit1)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 0.25, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA), main = "Approximated model")
scatter3D(x, y, z, pch = 18, cex = 0.01, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit1 = fit1points), main = "Approximated model")
scatter3D(x, y, z, pch = 18, cex = 0.01, 
          theta = 270, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit1 = fit1points), main = "approximated model")
summary(fit1)                 #log(dipo) seems to make some difference, interaction term can be excluded (high p-value), Dipo has a very small coefficient
par(mfrow=c(1,2))
plot(LogY,fit$residuals)
plot(LogY,fit1$residuals)      #residuals distribution appears slighlty better
#still problems at high values of logY
detach(DATA_0.50_NO)
DATA<-DATA[-which(is.na(DATA$UsableVS30)),]
attach(DATA_0.50_NO)detach(DATA_0.50_NO)

#TRying a new model: introduction of Mh=5.1
mh=4
x<-ML
y <- Dipo
k <- UsableVS30
z <- LogY
fit2 <- lm(z ~ ifelse((x-mh)<=0,x-mh,0)+ifelse((x-mh)>0,x-mh,0) + y +log10(y)+log10(k/800))
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
k.pred <- seq(min(k), max(k), length.out = grid.lines)
xyk <- expand.grid( x = x.pred, y = y.pred, k=k.pred)
z.pred <- matrix(predict(fit2, newdata = xyk), 
                 nrow = grid.lines, ncol = grid.lines)
# fit2ted points for droplines to surface
fit2points <- predict(fit2)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 0.25, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA), main = "RegPlane")
scatter3D(x, y, z, pch = 18, cex = 0.01, 
          theta = -60, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit2 = fit2points), main = "RegPlane")
scatter3D(x, y, z, pch = 18, cex = 0.01, 
          theta = 0, phi = 20, ticktype = "detailed",
          xlab = "ML", ylab = "Dipo", zlab = "LogY",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit2 = fit2points), main = "RegPlane")
summary(fit2)                 #log(dipo) seems to make some difference, interaction term can be excluded (high p-value), Dipo has a very small coefficient
par(mfrow=c(1,3))
plot(LogY,fit$residuals)
plot(LogY,fit1$residuals)
plot(LogY,fit2$residuals)   #actually the fitting it's not so much different
sum(abs(fit$residuals))
sum(abs(fit1$residuals))
sum(abs(fit2$residuals))





################ WRITE THE NEW TABLE #################

write.table(DATA,file=" FAS_000_50_NoOutlier.csv",append=FALSE, quote=FALSE,sep=",", eol="\n",na="NA",
            dec = ".", row.names = FALSE, col.names = TRUE)
