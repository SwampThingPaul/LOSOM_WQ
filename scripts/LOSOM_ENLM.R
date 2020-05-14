## 
## LOSOM Water Quality Subteam
## Estuary Nutrient Load Model
##
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
library(zyp)
library(mblm)
library(rkt)
library(car)
library(pedometrics)
library(MASS)
library(relaimpo)
library(lmtest)

#Paths
wd="C:/Julian_LaCie/_Github/LOSOM_WQ"

paths=paste0(wd,c("/Plots/","/Export/","/Data/"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

##
dates=as.Date(c("1998-05-01","2019-04-30"))
# CRE Hydro ---------------------------------------------------------------
# Discharge
q.dbkeys=data.frame(SITE=c("S79","S78","S77"),DBKEY=c("00865","DJ236","15635"))
q.cre.dat=data.frame()
for(i in 1:nrow(q.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],q.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(q.dbkeys$DBKEY[i])
  q.cre.dat=rbind(q.cre.dat,tmp)
  print(i)
}
q.cre.dat=merge(q.cre.dat,q.dbkeys,"DBKEY")
q.cre.dat$Date.EST=date.fun(q.cre.dat$Date)

q.cre.dat.xtab=cast(q.cre.dat,Date.EST~SITE,value="Data.Value",mean)
q.cre.dat.xtab$month=format(q.cre.dat.xtab$Date,"%m")
q.cre.dat.xtab$CY=format(q.cre.dat.xtab$Date,"%Y")
q.cre.dat.xtab$monCY=with(q.cre.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
q.cre.dat.xtab$WY=WY(q.cre.dat.xtab$Date.EST)

q.cre.dat.xtab$S77=with(q.cre.dat.xtab,ifelse(S77<0,0,S77))
q.cre.dat.xtab$S78=with(q.cre.dat.xtab,ifelse(S78<0,0,S78))
q.cre.dat.xtab$S78QFreq=as.numeric(q.cre.dat.xtab$S78>0)
q.cre.dat.xtab$C43.basin.in=with(q.cre.dat.xtab,ifelse(S79<S77,0,S79-S77))
q.cre.dat.xtab$basin.ratio=with(q.cre.dat.xtab,ifelse(S79==0,NA,C43.basin.in/S79))

#From Doering Model
q.cre.dat.xtab$LakeQ=with(q.cre.dat.xtab,ifelse(S79-S77>0,S77,S79))
# sanity check
# plot(LakeQ~Date.EST,q.cre.dat.xtab)
# subset(q.cre.dat.xtab,Date.EST==date.fun("2009-05-04"))

# Monthly aggregation

# Stage data
stg.dbkeys=data.frame(SITE=c("S79_H","S79_H","S235_T","S235_T","S77_H","S77_H"),DBKEY=c("00864","AN786","15566","38259","J8188","00852"))
stg.dbkeys=subset(stg.dbkeys,DBKEY!="00864");#not sure of the datum (NGVD vs NAVD)
# S235T in this analysis is a surragate for S77T
# S77H was added at the recommendation of SFWMD - explore the use of "Lake" Stage.
stg.cre.dat=data.frame()
for(i in 1:nrow(stg.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],stg.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(stg.dbkeys$DBKEY[i])
  stg.cre.dat=rbind(stg.cre.dat,tmp)
  print(i)
}
stg.cre.dat=merge(stg.cre.dat,stg.dbkeys,"DBKEY")
stg.dat.da.xtab=cast(stg.cre.dat,Date~SITE,value="Data.Value",mean)
stg.dat.da.xtab$WY=WY(stg.dat.da.xtab$Date)
stg.dat.da.xtab$month=format(stg.dat.da.xtab$Date,"%m")
stg.dat.da.xtab$CY=format(stg.dat.da.xtab$Date,"%Y")
stg.dat.da.xtab$monCY=with(stg.dat.da.xtab,date.fun(paste(CY,month,"01",sep="-")))
stg.dat.da.xtab$grad=with(stg.dat.da.xtab,(S235_T-S79_H))
#stg.dat.da.xtab$grad2=with(stg.dat.da.xtab,(S77_H-S79_H)); # Limited daily stage POR
#stg.dat.da.xtab$grad=with(stg.dat.da.xtab,(S235_T-S79_H)/213356.3)

# Sanity Check
# head(stg.dat.da.xtab)
# plot(grad2~Date,stg.dat.da.xtab,type="l")
# with(stg.dat.da.xtab,lines(Date,grad,col="red"))

# Monthly aggregation

# CRE WQ ------------------------------------------------------------------
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179),param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla"))
wq.sites=c("S79","S77")
wq.dat=DBHYDRO_WQ(dates[1],dates[2],wq.sites,params$Test.Number)
wq.dat=merge(wq.dat,params,"Test.Number")
unique(wq.dat$Collection.Method)
wq.dat=subset(wq.dat,Collection.Method=="G")

# POR Sampling plot
plot.order=data.frame(param=c("NOx","NH4","TKN","TN","SRP","TP","Chla"),
                      param.plt.order=c(1,2,3,4,5,6,7),
                      param.lab=c("NO\u2093","NH\u2084","TKN","TN","SRP","TP","Chl-a"))
wq.dat=merge(wq.dat,plot.order,"param")
unique(wq.dat$param)

ylim.val=c(1,nrow(plot.order));ymaj=seq(ylim.val[1],ylim.val[2],1)
xlim.val=date.fun(dates);xmaj=seq(xlim.val[1],xlim.val[2],"5 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
#tiff(filename=paste0(plot.path,"S79_wq_monitoring.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_wq_monitoring.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.5),oma=c(2,1,0.5,0.5));
tmp=subset(wq.dat,Station.ID==wq.sites[1])

plot(xlim.val,0:1,type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA)
abline(v=xmin,lty=2,col="grey")
with(hfr.clim,sapply(seq_along(prec),function(k,z) mean(z[1:k],na.rm=T),z=prec))
for(i in 1:length(unique(wq.dat$param))){
  with(subset(tmp,param==unique(param)[i]),points(Date.EST,param.plt.order,pch="|",col=adjustcolor("red",0.5)))
}
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymaj,plot.order$param.lab);box(lwd=1)
mtext(side=3,wq.sites[1])
mtext(side=1,line=1.75,"Date (Month-Year)")
mtext(side=2,line=3,"Parameter")
dev.off()

# WQ crosstab
wq.dat.xtab=cast(wq.dat,Station.ID+Date.EST~param,value="HalfMDL",mean)
wq.dat.xtab$DIN=with(wq.dat.xtab,NH4+NOx)
wq.dat.xtab$TN=with(wq.dat.xtab, TN_Combine(NOx,TKN,TN))
wq.dat.xtab$WY=WY(wq.dat.xtab$Date.EST)
wq.dat.xtab$month=format(wq.dat.xtab$Date.EST,"%m")
wq.dat.xtab$CY=format(wq.dat.xtab$Date.EST,"%Y")
wq.dat.xtab$monCY=with(wq.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))

# Reversal Evaluation
wq.dat.xtab$TPReversal=with(wq.dat.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
wq.dat.xtab$TNReversal=with(wq.dat.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

sum(wq.dat.xtab$TNReversal,na.rm=T)
subset(wq.dat.xtab,TNReversal==T)
sum(wq.dat.xtab$TPReversal,na.rm=T)

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,wq.dat.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,wq.dat.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

dev.off()

# Monthly aggregation

# Doering Model
WYs=seq(2010,2014,1)
#WYs=seq(1999,2019,1)
wq.vars=c("Station.ID","WY", "Date.EST", "Chla", "NH4", "NOx", "SRP", "TKN","TN", "TP", "DIN")
doering.dat1=merge(subset(q.cre.dat.xtab,WY%in%WYs),subset(wq.dat.xtab,WY%in%WYs&Station.ID=="S79")[,wq.vars],c("Date.EST","WY"))
head(doering.dat1)

doering.dat1$TP_kgd=with(doering.dat1,Load.Calc.kg(S79,TP))
doering.dat1$TN_kgd=with(doering.dat1,Load.Calc.kg(S79,TN))

## TP Model
doering.mod.TP=lm(TP_kgd~C43.basin.in+LakeQ,doering.dat1)

summary(doering.mod.TP)
gvlma::gvlma(doering.mod.TP)
layout(matrix(1:4,2,2));plot(doering.mod.TP)
shapiro.test(residuals(doering.mod.TP));hist(residuals(doering.mod.TP))

layout(matrix(1:2,1,2))
#Fit-Mean
fm.spread=ecdf_fun(doering.mod.TP$fitted.values-mean(doering.mod.TP$fitted.values))
plot(value~proportion,fm.spread);mtext(side=3,"Fit - Mean")
#Residual
rs.spread=ecdf_fun(doering.mod.TP$residuals)
plot(value~proportion,rs.spread);mtext(side=3,"Residual")

## TN Model
doering.mod.TN=lm(TN_kgd~C43.basin.in+LakeQ,doering.dat1)
summary(doering.mod.TN)
gvlma::gvlma(doering.mod.TN)
layout(matrix(1:4,2,2));plot(doering.mod.TN)
shapiro.test(residuals(doering.mod.TN));hist(residuals(doering.mod.TN))

layout(matrix(1:2,1,2))
#Fit-Mean
fm.spread=ecdf_fun(doering.mod.TN$fitted.values-mean(doering.mod.TN$fitted.values))
plot(value~proportion,fm.spread);mtext(side=3,"Fit - Mean")
#Residual
rs.spread=ecdf_fun(doering.mod.TN$residuals)
plot(value~proportion,rs.spread);mtext(side=3,"Residual")
