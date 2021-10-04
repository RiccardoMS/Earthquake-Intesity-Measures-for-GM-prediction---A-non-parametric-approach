#### thin plate spline visualization
LogY<-DATA_NP$LogY
Dipo<-DATA_NP$Dipo
ML<-DATA_NP$ML
model_gam_fic=gam(LogY ~ s(Dipo,ML, bs='tp'))
summary(model_gam_fic)

grid.dipo<-seq(range(DATA_NP$Dipo)[1],range(DATA_NP$Dipo)[2], length.out = 100)
grid.ml<-seq(range(DATA_NP$ML)[1],range(DATA_NP$ML)[2], length.out = 100)
grid<-expand.grid(grid.dipo,grid.ml)
names(grid)<-c("Dipo","ML")
newdata=data.frame(Dipo=grid$Dipo, ML=grid$ML)
pred=predict(model_gam_fic,newdata=newdata,type='response')
library(rgl)
persp3d(grid.dipo,grid.ml,pred,col='dodgerblue2',xlab = 'Dipo',ylab='ML', zlab='Thin-plate')
points3d(DATA_NP$Dipo,DATA_NP$ML,DATA_NP$LogY,col='darkgrey',size=1.5)
