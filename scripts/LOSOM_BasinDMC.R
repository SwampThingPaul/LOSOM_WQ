## 
## LOSOM Water Quality Subteam
## Estuary Nutrient Load Model
##   Estuary-Basin Double Mass Curve
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()


#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)

# Additional Analysis libraries
library(Hmisc)

library(zyp)
library(mblm)
library(rkt)
library(car)
library(pedometrics)
library(MASS)
library(relaimpo)
library(lmtest)

library(kableExtra)
#Paths
wd="C:/Julian_LaCie/_Github/LOSOM_WQ"

paths=paste0(wd,c("/Plots/","/Export/","/Data/"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

##
dates=as.Date(c("1978-05-01","2019-04-30"))

# Rainfall Data -----------------------------------------------------------
# STL/C44 Basin
# S308R=DBHYDRO_daily(dates[1],dates[2],c("06239","06119","16588"))
# plot(Data.Value~Date,S308R)
# range(S308R$Data.Value,na.rm=T)
# range(S308R$Date,na.rm=T)
# ddply(S308R,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))
# 
# S135R=DBHYDRO_daily(dates[1],dates[2],c("05849","K8637"))
# plot(Data.Value~Date,S135R)
# range(S135R$Data.Value,na.rm=T)
# range(S135R$Date,na.rm=T)
# ddply(S135R,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))
# 
# S80R=DBHYDRO_daily(dates[1],dates[2],c("06075","06237","16618"))
# plot(Data.Value~Date,S80R)
# range(S80R$Data.Value,na.rm=T)
# range(S80R$Date,na.rm=T)
# ddply(S80R,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))
# 
# PRATTR=DBHYDRO_daily(dates[1],dates[2],"06122")
# plot(Data.Value~Date,PRATTR)
# range(PRATTR$Data.Value,na.rm=T)
# range(PRATTR$Date,na.rm=T)
# 
# ACRA=DBHYDRO_daily(dates[1],dates[2],c("SX445","VM862"))
# plot(Data.Value~Date,ACRA)
# range(ACRA$Data.Value,na.rm=T)
# range(ACRA$Date,na.rm=T)
# ddply(ACRA,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))

C44Basin_rainfall.sites=data.frame(DBKEY=c(c("06239","06119","16588"),c("05849","K8637"),c("06075","06237","16618"),"06122",c("SX445","VM862")),
                                   SITE=c(rep("S308R",3),rep("S135R",2),rep("S80R",3),"PRATT",rep("ACRA2",2)))
C44Basin_rainfall.sites$BASIN="C44"

# CRE/C43 Basin
# PALMDALE=DBHYDRO_daily(dates[1],dates[2],c("06093","15786"))
# plot(Data.Value~Date,PALMDALE)
# range(PALMDALE$Data.Value,na.rm=T)
# range(PALMDALE$Date,na.rm=T)
# ddply(PALMDALE,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))
# 
# WHIDDEN=DBHYDRO_daily(dates[1],dates[2],c("VN424","15465"))
# plot(Data.Value~Date,WHIDDEN)
# range(WHIDDEN$Data.Value,na.rm=T)
# range(WHIDDEN$Date,na.rm=T)
# ddply(WHIDDEN,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))
# 
# S78R=DBHYDRO_daily(dates[1],dates[2],c("06221","16625","06243","15495"))
# plot(Data.Value~Date,S78R)
# range(S78R$Data.Value,na.rm=T)
# range(S78R$Date,na.rm=T)
# ddply(S78R,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))
# 
# KERI_TOW=DBHYDRO_daily(dates[1],dates[2],c("06083"))
# plot(Data.Value~Date,KERI_TOW)
# range(KERI_TOW$Data.Value,na.rm=T)
# range(KERI_TOW$Date,na.rm=T)
# 
# DEVILS=DBHYDRO_daily(dates[1],dates[2],c("05953","06079","64082"))
# plot(Data.Value~Date,DEVILS)
# range(DEVILS$Data.Value,na.rm=T)
# range(DEVILS$Date,na.rm=T)
# ddply(DEVILS,"DBKEY",summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))

C43Basin_rainfall.sites=data.frame(DBKEY=c(c("06093","15786"),c("VN424","15465"),c("06221","16625","06243","15495"),"06083",c("05953","06079","64082")),
                                   SITE=c(rep("PALMDALE",2),rep("WHIDDEN",2),rep("S78R",4),"KERI_TOW",rep("DEVILS",3)))
C43Basin_rainfall.sites$BASIN="C43"

rainfall.sites=rbind(C44Basin_rainfall.sites,C43Basin_rainfall.sites)
q.dbkeys=data.frame(SITE=c("S79","S80"),DBKEY=c("00865","JW224"),BASIN=c("C43","C44"))


# -------------------------------------------------------------------------
rf.dat=data.frame()
pb=txtProgressBar(1,nrow(rainfall.sites),1,style=3)

for(i in 1:nrow(rainfall.sites)){
  tmp=DBHYDRO_daily(dates[1],dates[2],rainfall.sites$DBKEY[i])
  tmp$DBKEY=as.character(rainfall.sites$DBKEY[i])
  rf.dat=rbind(tmp,rf.dat)
  setTxtProgressBar(pb, i)
}
rf.dat=merge(rf.dat,rainfall.sites,"DBKEY")
rf.dat$Date.EST=date.fun(rf.dat$Date)
rf.dat.da=ddply(rf.dat,c("SITE","BASIN","Date.EST"),summarise,mean.RF=mean(in.to.cm(Data.Value),na.rm=T))
rf.dat.da2=ddply(rf.dat.da,c("BASIN","Date.EST"),summarise,mean.RF=mean(mean.RF,na.rm=T),N.val=N(mean.RF))
rf.dat.da2$WY=WY(rf.dat.da2$Date.EST)

ddply(rf.dat,c("SITE","BASIN"),summarise,min.dat=min(Date.EST),max.dat=max(Date.EST))

rf.dat.WY=ddply(rf.dat.da2,c("BASIN","WY"),summarise,TRF.cm=sum(mean.RF,na.rm=T))
rf.dat.WY$cum.rf=with(rf.dat.WY,ave(TRF.cm,BASIN,FUN=function(x)cumsum(x)))

plot(TRF.cm~WY,subset(rf.dat.WY,BASIN=="C44"),type="l")
with(subset(rf.dat.WY,BASIN=="C43"),lines(WY,TRF.cm,col="red"))


q.dat=data.frame()
pb=txtProgressBar(1,nrow(q.dbkeys),1,style=3)
for(i in 1:nrow(q.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],q.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(q.dbkeys$DBKEY[i])
  q.dat=rbind(tmp,q.dat)
  setTxtProgressBar(pb, i)
}
q.dat=merge(q.dat,q.dbkeys,"DBKEY")

q.dat$Date.EST=date.fun(q.dat$Date)
q.dat$WY=WY(q.dat$Date.EST)

q.dat.WY=ddply(q.dat,c("BASIN","WY"),summarise,TFlow=sum(cfs.to.km3d(Data.Value),na.rm=T))
q.dat.WY$cum.q=with(q.dat.WY,ave(TFlow,BASIN,FUN=function(x)cumsum(x)))

q.rf=merge(q.dat.WY,rf.dat.WY,c("BASIN","WY"))

plot(cum.q~cum.rf,subset(q.rf,BASIN=="C43"))
plot(cum.q~cum.rf,subset(q.rf,BASIN=="C44"))

library(segmented)

c43.rq=lm(cum.q~cum.rf,subset(q.rf,BASIN=="C43"))
layout(matrix(1:4,2,2));plot(c43.rq)
dev.off()
shapiro.test(residuals(c43.rq))

c43.rq.seg=segmented(c43.rq,seg.Z=~cum.rf,psi=c(1150,2000,3900,4900))
summary(c43.rq.seg)
plot(cum.q~cum.rf,subset(q.rf,BASIN=="C43"))
plot(c43.rq.seg,add=T)

resvar=sum(c43.rq.seg$residuals^2)/c43.rq.seg$df.residual
mss= if (attr(c43.rq.seg$terms, "intercept")){sum((c43.rq.seg$fitted.values - mean(c43.rq.seg$fitted.values))^2)}else{sum(c43.rq.seg$fitted.values^2)}
p <- c43.rq.seg$rank
df.int<-if (attr(c43.rq.seg$terms, "intercept")){1L}else{0L}
Fval=(mss/(p - df.int))/resvar;#F-value
numdf=p - df.int;#numdf
denf=c43.rq.seg$df.residual;#denf
pval=pf((mss/(p - df.int))/resvar,p - df.int,c43.rq.seg$df.residual,lower.tail = F);#pvalue


##
c44.rq=lm(cum.q~cum.rf,subset(q.rf,BASIN=="C44"))
layout(matrix(1:4,2,2));plot(c44.rq)
dev.off()
shapiro.test(residuals(c44.rq))

c44.rq.seg=segmented(c44.rq,seg.Z=~cum.rf,psi=c(2500,3200,4200,5050))
summary(c44.rq.seg)
plot(cum.q~cum.rf,subset(q.rf,BASIN=="C44"))
plot(c44.rq.seg,add=T)

BK.pt=c44.rq.seg$psi[,2]
BK.pt.SE=c44.rq.seg$psi[,3]
RSE=summary(c44.rq.seg)$sigma
R2=summary(c44.rq.seg)$r.squared
resvar=sum(c44.rq.seg$residuals^2)/c44.rq.seg$df.residual
mss= if (attr(c44.rq.seg$terms, "intercept")){sum((c44.rq.seg$fitted.values - mean(c44.rq.seg$fitted.values))^2)}else{sum(c44.rq.seg$fitted.values^2)}
p <- c44.rq.seg$rank
df.int<-if (attr(c44.rq.seg$terms, "intercept")){1L}else{0L}
Fval=(mss/(p - df.int))/resvar;#F-value
numdf=p - df.int;#numdf
denf=c44.rq.seg$df.residual;#denf
pval=pf((mss/(p - df.int))/resvar,p - df.int,c44.rq.seg$df.residual,lower.tail = F);#pvalue


#tiff(filename=paste0(plot.path,"C43_DMC.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));

xlim.val=c(0,6000);by.x=1000;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,80);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(cum.q~cum.rf,q.rf,type="n",ylab=NA,xlab=NA,axes=F,xlim=xlim.val,ylim=ylim.val)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(q.rf,BASIN=="C43"),pt_line(cum.rf,cum.q,2,"dodgerblue1",1,21,"dodgerblue1",pt.lwd=0.1))
x.val=with(subset(q.rf,BASIN=="C43"),seq(min(cum.rf),max(cum.rf),length.out=50))
seg.mod.pred=predict(c43.rq.seg,data.frame(cum.rf=x.val),interval="confidence")
lines(x.val,seg.mod.pred[,1],col=adjustcolor("indianred1",0.5),lwd=3)
abline(v=c43.rq.seg$psi[,2],lty=2)
yrs=data.frame()
for(i in 1:length(c43.rq.seg$psi[,2])){
  yrs=rbind(yrs,data.frame(WY=max(subset(q.rf,BASIN=="C43"&cum.rf<c43.rq.seg$psi[i,2])$WY)))
}
text(c43.rq.seg$psi[,2],rep(ylim.val[2],length(yrs)),yrs$WY,cex=0.75)  
axis_fun(1,line=-0.5,xmaj,xmin,xmaj*0.01)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,"C-43 Basin")
mtext(side=2,line=2.25,"Cumulative Discharge (km\u00B3 Yr\u207B\u00b9)")
mtext(side=1,line=1.75,"Cumulative Rainfall (m Yr\u207B\u00b9)")
dev.off()

#tiff(filename=paste0(plot.path,"C43C44_DMC.tiff"),width=6.5,height=3.25,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
layout(matrix(1:2,1,2))

xlim.val=c(0,6000);by.x=1000;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,80);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(cum.q~cum.rf,q.rf,type="n",ylab=NA,xlab=NA,axes=F,xlim=xlim.val,ylim=ylim.val)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(q.rf,BASIN=="C43"),pt_line(cum.rf,cum.q,2,"dodgerblue1",1,21,"dodgerblue1",pt.lwd=0.1))
x.val=with(subset(q.rf,BASIN=="C43"),seq(min(cum.rf),max(cum.rf),length.out=50))
seg.mod.pred=predict(c43.rq.seg,data.frame(cum.rf=x.val),interval="confidence")
lines(x.val,seg.mod.pred[,1],col=adjustcolor("indianred1",0.5),lwd=3)
abline(v=c43.rq.seg$psi[,2],lty=2)
yrs=data.frame()
for(i in 1:length(c43.rq.seg$psi[,2])){
yrs=rbind(yrs,data.frame(WY=max(subset(q.rf,BASIN=="C43"&cum.rf<c43.rq.seg$psi[i,2])$WY)))
}
text(c43.rq.seg$psi[,2],rep(ylim.val[2],length(yrs)),yrs$WY,cex=0.75)  
axis_fun(1,line=-0.5,xmaj,xmin,xmaj*0.01)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,"C-43 Basin")
mtext(side=2,line=2.5,"Cumulative Discharge (km\u00B3 Yr\u207B\u00b9)")

ylim.val=c(0,20);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(cum.q~cum.rf,q.rf,type="n",ylab=NA,xlab=NA,axes=F,xlim=xlim.val,ylim=ylim.val)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(q.rf,BASIN=="C44"),pt_line(cum.rf,cum.q,2,"dodgerblue1",1,21,"dodgerblue1",pt.lwd=0.1))
x.val=with(subset(q.rf,BASIN=="C44"),seq(min(cum.rf),max(cum.rf),length.out=50))
seg.mod.pred=predict(c44.rq.seg,data.frame(cum.rf=x.val),interval="confidence")
lines(x.val,seg.mod.pred[,1],col=adjustcolor("indianred1",0.5),lwd=3)
abline(v=c44.rq.seg$psi[,2],lty=2)
yrs=data.frame()
for(i in 1:length(c44.rq.seg$psi[,2])){
  yrs=rbind(yrs,data.frame(WY=max(subset(q.rf,BASIN=="C44"&cum.rf<c44.rq.seg$psi[i,2])$WY)))
}
text(c44.rq.seg$psi[,2],rep(ylim.val[2],length(yrs)),yrs$WY,cex=0.75)  
axis_fun(1,line=-0.5,xmaj,xmin,xmaj*0.01)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,"C-44 Basin")
mtext(side=1,line=0.5,"Cumulative Rainfall (m Yr\u207B\u00b9)",outer=T)
dev.off()
