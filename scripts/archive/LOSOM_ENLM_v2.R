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

library(dataRetrieval)
library(RCurl)

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

# q.cre.dat.xtab$S77.30d=with(q.cre.dat.xtab,c(rep(NA,29),rollsum(S77,k=30,na.rm=T)))
# q.cre.dat.xtab$S78.30d=with(q.cre.dat.xtab,c(rep(NA,29),rollsum(S78,k=30,na.rm=T)))
# q.cre.dat.xtab$S79.30d=with(q.cre.dat.xtab,c(rep(NA,29),rollsum(S79,k=30,na.rm=T)))
# q.cre.dat.xtab$C43basin.30d=with(q.cre.dat.xtab,ifelse(S79.30d<S77.30d,0,S79.30d-S77.30d))

##
#From Doering Model
q.cre.dat.xtab$LakeQ=with(q.cre.dat.xtab,ifelse(S79-S77>0,S77,S79))
# sanity check
# plot(LakeQ~Date.EST,q.cre.dat.xtab)
# subset(q.cre.dat.xtab,Date.EST==date.fun("2009-05-04"))

# Monthly aggregation
q.cre.dat.xtab.mon=ddply(q.cre.dat.xtab,c("monCY","WY"),summarise,S77=sum(cfs.to.km3d(S77),na.rm=T),S78=sum(cfs.to.km3d(S78),na.rm=T),S79=sum(cfs.to.km3d(S79),na.rm=T),S78Qfreq=sum(S78QFreq,na.rm=T))
q.cre.dat.xtab.mon$C43Basin=with(q.cre.dat.xtab.mon,ifelse(S79<S77,0,S79-S77))
q.cre.dat.xtab.mon$LakeQ=with(q.cre.dat.xtab.mon,ifelse(S79-S77>0,S77,S79))
q.cre.dat.xtab.mon$basin.q.ratio=with(q.cre.dat.xtab.mon,C43Basin/S79)

# Stage data
lakeo.stg.dbkey=data.frame(SITE=c("LakeO","LakeO"),DBKEY=c("15611","06832"),priority=c("P1","P2"));#15611 Pref DBKEY
lakeO.stg=DBHYDRO_daily(dates[1],dates[2],lakeo.stg.dbkey$DBKEY[1])
plot(Data.Value~Date,subset(lakeO.stg,DBKEY=="15611"))

lakeO.stg$Date.EST=date.fun(lakeO.stg$Date)

lakeO.stg$WY=WY(lakeO.stg$Date)
lakeO.stg$month=format(lakeO.stg$Date,"%m")
lakeO.stg$CY=format(lakeO.stg$Date,"%Y")
lakeO.stg$monCY=with(lakeO.stg,date.fun(paste(CY,month,"01",sep="-")))

# Monthly aggregation
lakeO.stg.mon=ddply(lakeO.stg,c("monCY","WY"),summarise,
                          LakeO=mean(Data.Value,na.rm=T))
plot(LakeO~monCY,lakeO.stg.mon)

# CRE WQ ------------------------------------------------------------------
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7),param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp"))
wq.sites=c("S79","S77")
wq.dat=DBHYDRO_WQ(dates[1],dates[2],wq.sites,params$Test.Number)
wq.dat=merge(wq.dat,params,"Test.Number")
unique(wq.dat$Collection.Method)
wq.dat=subset(wq.dat,Collection.Method=="G")

# POR Sampling plot
plot.order=data.frame(param=c("NOx","NH4","TKN","TN","SRP","TP","Chla","Temp"),
                      param.plt.order=c(1,2,3,4,5,6,7,8),
                      param.lab=c("NO\u2093","NH\u2084","TKN","TN","SRP","TP","Chl-a","Temp"))
wq.dat=merge(wq.dat,plot.order,"param",all.x=T)
unique(wq.dat$param)

ylim.val=c(1,nrow(plot.order));ymaj=seq(ylim.val[1],ylim.val[2],1)
xlim.val=date.fun(dates);xmaj=seq(xlim.val[1],xlim.val[2],"5 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
#tiff(filename=paste0(plot.path,"S79_wq_monitoring.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S79_wq_monitoring.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.5),oma=c(2,1,0.5,0.5));
tmp=subset(wq.dat,Station.ID==wq.sites[1])

plot(xlim.val,0:1,type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA)
abline(v=xmin,lty=2,col="grey")
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
wq.dat.xtab$TON=with(wq.dat.xtab,TN-DIN)
wq.dat.xtab$WY=WY(wq.dat.xtab$Date.EST)
wq.dat.xtab$month=format(wq.dat.xtab$Date.EST,"%m")
wq.dat.xtab$CY=format(wq.dat.xtab$Date.EST,"%Y")
wq.dat.xtab$monCY=with(wq.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
wq.dat.xtab$TNTP=with(wq.dat.xtab, (TN/14.00)/(TP/30.97))

boxplot(TON~Station.ID,wq.dat.xtab)

unique(wq.dat.xtab$Station.ID)
unique(q.cre.dat$SITE)
q.cre.dat$Q.cfs=q.cre.dat$Data.Value

subset(wq.dat.xtab,Station.ID=="S79"&month=="06"&CY==2011)
subset(q.cre.dat.xtab,Date.EST%in%date.fun(c("2011-06-01","2011-06-08")))

wq.dat.xtab=merge(wq.dat.xtab,q.cre.dat[,c("SITE","Date.EST","Q.cfs")],by.x=c("Station.ID","Date.EST"),by.y=c("SITE","Date.EST"),all.x=T)
range(wq.dat.xtab$Q.cfs,na.rm=T)

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
wq.dat.xtab.mon=ddply(subset(wq.dat.xtab,Station.ID=="S79"&abs(Q.cfs)>0),c("WY","monCY"),summarise,mean.TP=mean(TP,na.rm=T),mean.TN=mean(TN,na.rm=T),mean.DIN=mean(DIN,na.rm=T),mean.TON=mean(TON,na.rm=T),mean.TNTP=mean(TNTP,na.rm=T))
wq.dat.xtab.mon
#

# WYs=seq(2010,2014,1)
# #WYs=seq(1999,2019,1)
# wq.vars=c("Station.ID","WY", "Date.EST", "Chla", "NH4", "NOx", "SRP", "TKN","TN", "TP", "DIN")
# doering.dat1=merge(subset(q.cre.dat.xtab,WY%in%WYs),subset(wq.dat.xtab,WY%in%WYs&Station.ID=="S79")[,wq.vars],c("Date.EST","WY"))
# head(doering.dat1)
# 
# doering.dat1$TP_kgd=with(doering.dat1,Load.Calc.kg(S79,TP))
# #doering.dat1$TP_kgd.detrend=doering.dat1$TP_kgd-mean(doering.dat1$TP_kgd,na.rm=T)
# doering.dat1$TN_kgd=with(doering.dat1,Load.Calc.kg(S79,TN))
# #
# ## TP Model
# set.seed(1)
# tr.index=sample(1:nrow(doering.dat1),nrow(doering.dat1)*0.7)
# 
# doering.mod.TP=lm(log(TP_kgd)~C43.basin.in+LakeQ,doering.dat1[tr.index,])
# summary(doering.mod.TP)
# gvlma::gvlma(doering.mod.TP)
# layout(matrix(1:4,2,2));plot(doering.mod.TP)
# shapiro.test(residuals(doering.mod.TP));hist(residuals(doering.mod.TP))
#
# layout(matrix(1:2,1,2))
# #Fit-Mean
# fm.spread=ecdf_fun(doering.mod.TP$fitted.values-mean(doering.mod.TP$fitted.values))
# plot(value~proportion,fm.spread);mtext(side=3,"Fit - Mean")
# #Residual
# rs.spread=ecdf_fun(doering.mod.TP$residuals)
# plot(value~proportion,rs.spread);mtext(side=3,"Residual")
#
# ## TN Model
# doering.mod.TN=lm(TN_kgd~C43.basin.in+LakeQ,doering.dat1)
# summary(doering.mod.TN)
# gvlma::gvlma(doering.mod.TN)
# layout(matrix(1:4,2,2));plot(doering.mod.TN)
# shapiro.test(residuals(doering.mod.TN));hist(residuals(doering.mod.TN))
#
# layout(matrix(1:2,1,2))
# #Fit-Mean
# fm.spread=ecdf_fun(doering.mod.TN$fitted.values-mean(doering.mod.TN$fitted.values))
# plot(value~proportion,fm.spread);mtext(side=3,"Fit - Mean")
# #Residual
# rs.spread=ecdf_fun(doering.mod.TN$residuals)
# plot(value~proportion,rs.spread);mtext(side=3,"Residual")



# Concentration Discharge Relationship ------------------------------------
mean(subset(q.cre.dat.xtab,WY>2009)$S79,na.rm=T)
mean(q.cre.dat.xtab$S79,na.rm=T)
q.cre.dat.xtab$S79.QQmean=with(q.cre.dat.xtab,S79/mean(subset(q.cre.dat.xtab,WY>2009)$S79,na.rm=T))

plot(S79~Date.EST,q.cre.dat.xtab,type="l")
abline(h=mean(subset(q.cre.dat.xtab,WY>2009)$S79,na.rm=T),lty=2,col="indianred1")
abline(h=mean(q.cre.dat.xtab$S79,na.rm=T),lty=2,col="indianred1")

wq.vars=c("Station.ID","WY", "Date.EST", "Chla", "NH4", "NOx", "SRP", "TKN","TN", "TP", "DIN","TON")
cq.dat=merge(subset(q.cre.dat.xtab,WY>2009),subset(wq.dat.xtab,WY>2009&Station.ID=="S79")[,wq.vars],c("Date.EST","WY"))
cq.dat$log.S79.QQmean=log(cq.dat$S79.QQmean)
cq.dat$log.S79=log(cq.dat$S79+1)
cq.dat$log.TP=log(cq.dat$TP)
cq.dat$log.TN=log(cq.dat$TN)
head(cq.dat)

plot(TP~S79.QQmean,cq.dat,log="xy")
with(cq.dat,cor.test(TP,S79.QQmean,method="spearman"))
cols=colorRampPalette(c("indianred1","forestgreen"))(10)
tmp=subset(cq.dat, is.infinite(log.S79.QQmean)==F&is.na(log.S79.QQmean)==F&is.na(log.TP)==F)
for(i in 1:10){
  TPcQ=mblm(log.TP~log.S79.QQmean,subset(tmp,WY==(2010:2019)[i]))
  #summary(TPcQ)
  x.val=with(subset(tmp,WY==(2010:2019)[i]),seq(min(log.S79.QQmean,na.rm=T),max(log.S79.QQmean,na.rm=T),length.out=50))
  TPcQ.pred=predict(TPcQ,data.frame(log.S79.QQmean=x.val),interval="confidence")
  lines(exp(x.val),exp(TPcQ.pred[,1]),col=cols[i],lwd=3)
  #lines(exp(x.val),exp(TPcQ.pred[,2]),lty=2)
  #lines(exp(x.val),exp(TPcQ.pred[,3]),lty=2)
}

plot(TN~S79.QQmean,cq.dat,log="xy")
with(cq.dat,cor.test(TN,S79.QQmean,method="spearman"))
cols=colorRampPalette(c("indianred1","forestgreen"))(10)
tmp=subset(cq.dat, is.infinite(log.S79.QQmean)==F&is.na(log.S79.QQmean)==F&is.na(log.TN)==F)
for(i in 1:10){
  TNcQ=mblm(log.TN~log.S79.QQmean,subset(tmp,WY==(2010:2019)[i]))
  #summary(TNcQ)
  x.val=with(subset(tmp,WY==(2010:2019)[i]),seq(min(log.S79.QQmean,na.rm=T),max(log.S79.QQmean,na.rm=T),length.out=50))
  TNcQ.pred=predict(TNcQ,data.frame(log.S79.QQmean=x.val),interval="confidence")
  lines(exp(x.val),exp(TNcQ.pred[,1]),col=cols[i],lwd=3)
  #lines(exp(x.val),exp(TNcQ.pred[,2]),lty=2)
  #lines(exp(x.val),exp(TNcQ.pred[,3]),lty=2)
}


with(cq.dat,cor.test(TN,S79,method="spearman"))
plot(TON~S79.QQmean,cq.dat,log="xy")
plot(DIN~S79.QQmean,cq.dat,log="xy")


WYs=2010:2019
cols=wesanderson::wes_palette("Zissou1",length(WYs),"continuous")
#tiff(filename=paste0(plot.path,"CRE_S79CQ.tiff"),width=6.5,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.1,0.1),oma=c(3,2.5,1,0.5));
layout(matrix(c(1,2,3,3),2,2),widths=c(1,0.25))
xlim.val=c(0.01,10);xmaj=log.scale.fun(xlim.val,"major");xmin=log.scale.fun(xlim.val,"minor")

ylim.val=c(0.05,0.500);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
plot(TP~S79.QQmean,cq.dat,log="xy",type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYs)){
  with(subset(cq.dat,WY==WYs[i]),points(S79.QQmean,TP,pch=21,lwd=0.01,col=adjustcolor("black",0.5),bg=adjustcolor(cols[i],0.5),cex=1.25))
}
axis_fun(1,xmaj,xmin,NA)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=2,line=2.5,"TP (\u03BCg L\u207B\u00B9)")
mtext(side=3,"S-79")

ylim.val=c(0.5,3);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
plot(TN~S79.QQmean,cq.dat,log="xy",type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYs)){
  with(subset(cq.dat,WY==WYs[i]),points(S79.QQmean,TN,pch=21,lwd=0.01,col=adjustcolor("black",0.5),bg=adjustcolor(cols[i],0.5),cex=1.25))
}
axis_fun(1,xmaj,xmin,format(xmaj))
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"TN (mg L\u207B\u00B9)")
mtext(side=1,line=2,expression(paste("Q"["S79"]," : Q"[bar("X")])))

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,WYs,
       pch=21,lty=NA,lwd=0.1,
       col=adjustcolor("black",0.5),
       pt.bg=adjustcolor(cols,0.5),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,title.adj = 0,title="Water Year")
dev.off()

# monthly conc model ------------------------------------------------------
cre.hydro.mon=merge(q.cre.dat.xtab.mon,lakeO.stg.mon,c("monCY","WY"))
cre.hydro.wq.mon=merge(cre.hydro.mon,wq.dat.xtab.mon,c("WY","monCY"))

cre.hydro.wq.mon$hydro.season=as.numeric(FL.Hydroseason(cre.hydro.wq.mon$monCY)=="B_Dry")
cre.hydro.wq.mon$month=as.numeric(format(cre.hydro.wq.mon$monCY,"%m"))
cre.hydro.wq.mon=merge(cre.hydro.wq.mon,data.frame(month=c(5:12,1:4),month.plot=1:12),"month",all.x=T)
cre.hydro.wq.mon=cre.hydro.wq.mon[order(cre.hydro.wq.mon$monCY),]
cre.hydro.wq.mon$mean.TP.ugL=cre.hydro.wq.mon$mean.TP*1000

cre.hydro.wq.mon=subset(cre.hydro.wq.mon,WY>2009)

vars=c("mean.TP","LakeO","C43Basin","S77","S78","basin.q.ratio","WY")
#tiff(filename=paste0(plot.path,"CRE_scatterplot_month.tiff"),width=7,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.1,0.1),oma=c(3,3.5,0.75,0.5));
layout(matrix(1:49,7,7))

params=c("mean.TP.ugL","mean.TN","S77","S79","S78Qfreq","C43Basin", "basin.q.ratio","LakeO")
summary(cre.hydro.wq.mon[,params])
axis.lab=c("S79 TP\n(\u03BCg L\u207B\u00B9)",
           "S79 TN\n(mg L\u207B\u00B9)",
           "Q S77\n(km\u00B3 mon\u207B\u00B9)",
           "Q S79\n(km\u00B3 mon\u207B\u00B9)",
           "Q S78>0 Freq\n(days)",
           "Q C43\n(km\u00B3 mon\u207B\u00B9)",
           "Basin Ratio\n(Ratio)",
           "Lake O Stg\n(Ft)")

for(j in 1:7){
  if(j!=1){for(k in 1:(j-1)){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}}
  
  params2=params[-1:-j]
  axis.lab2=axis.lab[-1:-j]
  lim.min=c(0,0,0,0,0,0,0,9.5)
  lim.max=c(300,3,0.50,1,40,0.60,1,17);by.val=c(150,1.5,0.25,0.5,15,0.3,0.5,2)
  for(i in 1:length(params2)){
    xlim.val=c(lim.min[j],lim.max[j]);by.x=by.val[j];xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
    ylim.val=c(lim.min[-1:-j][i],lim.max[-1:-j][i]);by.y=by.val[-1:-j][i];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
    xmaj.lab=if(max(xmaj)>1e4){xmaj/1000}else{xmaj}
    ymaj.lab=if(max(ymaj)>1e4){ymaj/1000}else{ymaj}
    plot(cre.hydro.wq.mon[,params[j]],cre.hydro.wq.mon[,params2[i]],xlim=xlim.val,ylim=ylim.val,axes=F,type="n",ylab=NA,xlab=NA)
    abline(h=ymaj,v=xmaj,lty=3,col="grey")
    points(cre.hydro.wq.mon[,params[j]],cre.hydro.wq.mon[,params2[i]],pch=21,bg=adjustcolor("dodgerblue1",0.25),col=adjustcolor("grey",0.5),lwd=0.01,cex=1)
    y.val=na.omit(cre.hydro.wq.mon[,c(params[j],params2[i])])[,params2[i]]
    x.val=na.omit(cre.hydro.wq.mon[,c(params[j],params2[i])])[,params[j]]
    mod=mblm(y.val~x.val)
    x.val.seq=seq(min(x.val),max(x.val),length.out=50)
    mod.pred=predict(mod,data.frame(x.val=x.val.seq),interval="confidence")
    lines(x.val.seq,mod.pred[,1],col="indianred1")
    lines(x.val.seq,mod.pred[,2],col="indianred1",lty=2)
    lines(x.val.seq,mod.pred[,3],col="indianred1",lty=2)
    if(i==length(params2)){axis_fun(1,xmaj,xmin,format(xmaj.lab),cex=0.8,line=-0.75)}else{axis_fun(1,xmaj,xmin,NA,cex=0.8,line=-0.5)}
    if(j==1){axis_fun(2,ymaj,ymin,format(ymaj.lab),cex=0.8)}else{axis_fun(2,ymaj,ymin,NA,cex=0.8)}
    box(lwd=1)
    if(j==1){mtext(side=2,line=2.25,cex=0.65,axis.lab2[i])}
  }
  mtext(side=1,line=2.5,cex=0.65,axis.lab[j])
}
dev.off()

cor.mon=rcorr(as.matrix(cre.hydro.wq.mon[,params]),type="spearman")
cor.mon$r=round(cor.mon$r,2)
cor.mon$r[upper.tri(cor.mon$r,diag=F)]=NA
cor.mon$P=with(cor.mon,ifelse(P<0.05,"<0.05",round(P,2)))
cor.mon$P[upper.tri(cor.mon$P,diag=F)]=NA
cor.mon2=matrix(with(cor.mon,paste0(r," (",P,")")),ncol=length(params))
diag(cor.mon2)=NA
cor.mon2[cor.mon2=="NA (NA)"]=NA
param.names=c("TP","TN","Q S77","Q S79","Q S78 Freq","Q C43","Basin Ratio","Lake O")
rownames(cor.mon2)=param.names
colnames(cor.mon2)=param.names

knitr::kable(cor.mon2[2:8,1:7],align=c("c"),escape=F)%>%
  kable_styling( full_width = F)

# Quick PCA
library(REdaS)
library(vegan)
#params=c("mean.TP.ugL","mean.TN","mean.DIN","mean.TON","S77","S78","S78Qfreq","C43Basin", "basin.q.ratio","LakeO")
params=c("mean.TP.ugL","mean.TN","S79","S77","S78","S78Qfreq","C43Basin", "basin.q.ratio","LakeO")

KMOS(na.omit(cre.hydro.wq.mon[,params]))
bart_spher(na.omit(cre.hydro.wq.mon[,params]))
my.pca=rda(na.omit(cre.hydro.wq.mon[,params]),scale=T)

plot(my.pca)

eig <- my.pca$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
eig.pca

#tiff(filename=paste0(plot.path,"CRE_monthly_PCA.tiff"),width=5,height=5.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,3,0.1,1),oma=c(0.5,2,0.75,0.5));
layout(matrix(c(1:3),3,1),heights=c(1,1,0.2))

ylim.val=c(0,5);by.y=1;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2)
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n")
abline(h=ymaj,lty=3,col="grey")
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n",add=T)
abline(h=1,lty=2,col="red",lwd=2)
axis_fun(1,line=-0.7,x,x,NA,0.7)
axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=2,line=2.25,"Eigenvalue")

ylim.val=c(0,110);by.y=25;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2);#set y limit and delineates the major and minor ticks
x=barplot(eig.pca$variance,ylim=ylim.val,col="white",border=0,yaxt="n")# inital plot to get the measurements
abline(h=ymaj,lty=3,col="grey")#makes vertical lines from y axis
x=barplot(eig.pca$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)# the real plot that matters
lines(x,eig.pca$cumvariance,col="indianred1",lwd=2)# adds the cumulative variance for each factor
points(x,eig.pca$cumvariance,pch=21,bg="indianred1",cex=2,lwd=0.01)
abline(h=80,lty=2,col="red",lwd=2)
axis_fun(1,x,x,seq(1,length(x),1),1)
axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=1,line=2,"Principal Components")
mtext(side=2,line=2.25,"Percentage of Variances")

plot(0:1,0:1,axes=F,ylab=NA,xlab=NA,type="n")
legend.text=c("Absolute","Cumulative");#helper vaiable for legend
pt.col=c("grey","indianred1")#helper vaiable for legend
legend(0.5,0,legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(0.5,0,legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1.55,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

scrs=scores(my.pca,display=c("sites","species"),choices=c(1,2,3));
rownames(scrs$species)
axis.lab=c("TP",
           "TN",
           "Q S79",
           "Q S77",
           "Q S78",
           "Q S78>0 Freq",
           "Q C43",
           "Basin Ratio",
           "Lake O Stg")

#tiff(filename=paste0(plot.path,"CRE_month_PCA_biplot.tiff"),width=4.5,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,3,0.1,0.5),oma=c(1.5,1.5,0.75,0.5));
#layout(matrix(1:2,2,1))

xlim.val=c(-1,2.5);by.x=1;xmaj=c(seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-2.5,2.5);by.y=1.25;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs$sites[,c(1,2)],pch=21,bg=adjustcolor("grey",0.20),col=adjustcolor("black",0.20),cex=1,lwd=0.1); #plots the points
arrows(0,0,scrs$species[,1],scrs$species[,2],length = 0.07, angle = 25, code = 2,col="indianred1",lwd=2);# makes the arrows
with(scrs,text(species[,1],species[,2],labels=axis.lab,cex=1,font=2,pos=4));#adds labels to the arrows; 
# with(scrs,text(species[2,1],species[2,2],labels=axis.lab[2],cex=1,font=2,pos=1));#adds labels to the arrows; 
# with(scrs,text(species[c(1,3,4,6,7,8),1],species[c(1,3,4,6,7,8),2],labels=axis.lab[c(1,3,4,6,7,8)],cex=1,font=2,pos=2));#adds labels to the arrows; 
# with(scrs,text(species[5,1],species[5,2],labels=axis.lab[5],cex=1,font=2,pos=1));#adds labels to the arrows; 
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.5,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2.5,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance

# plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
# abline(h=0,v=0,lty=3,col="grey");
# points(scrs$sites[,c(1,3)],pch=21,bg=adjustcolor("grey",0.50),col=adjustcolor("black",0.50),cex=1.5,lwd=0.1); #plots the points
# arrows(0,0,scrs$species[,1],scrs$species[,3],length = 0.07, angle = 25, code = 2,col="indianred1",lwd=2);# makes the arrows
# with(scrs,text(species[,1],species[,3],labels=axis.lab,cex=1,font=2,pos=2));#adds labels to the arrows;
# axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
# axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
# mtext(side=1,line=1.5,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
# mtext(side=2,line=2.5,paste0("PCA 3 (",round(eig.pca$variance[3],1),"%)"));#adds y axis label with percent variance
dev.off()


# Hydrodynamic WQ model ---------------------------------------------------

# TP Model ----------------------------------------------------------------

null.mod=lm(mean.TP~1,cre.hydro.wq.mon)
summary(null.mod)
layout(matrix(1:4,2,2));plot(null.mod)
shapiro.test(residuals(null.mod))
acf(residuals(null.mod))
#car::vif(null.mod)

#vars=c("mean.TP","grad","C43Basin","S78","S79","S77","basin.q.ratio","hydro.season")
#vars=c("mean.TP","grad","S79_H","S235_T","C43Basin","S77","S78","basin.q.ratio","S78Qfreq")
vars=c("mean.TP","C43Basin","S77","S78","basin.q.ratio","S78Qfreq","LakeO")
S79.TP.mod.all=lm(log(mean.TP)~.,na.omit(cre.hydro.wq.mon[,vars]))
#AIC model
S79.TP.mod.sw=stepAIC(S79.TP.mod.all,direction="both",trace=F)
S79.TP.mod.sw$anova
summary(S79.TP.mod.sw)
formula(S79.TP.mod.sw)
vif(S79.TP.mod.sw)

plot(mean.TP~monCY,cre.hydro.wq.mon,type="l")
with(cre.hydro.wq.mon,lines(monCY,lag(mean.TP,10),col="blue"))

S79.TP.mod=lm(log(mean.TP)~C43Basin+S77+LakeO,cre.hydro.wq.mon)
layout(matrix(1:4,2,2));plot(S79.TP.mod)
shapiro.test(residuals(S79.TP.mod));hist(residuals(S79.TP.mod))
acf(S79.TP.mod$residuals)
vif(S79.TP.mod)
summary(S79.TP.mod)

dev.off()
## Train vs test
# http://r-statistics.co/Linear-Regression.html
# using random 80% of data
# range(cre.hydro.wq.mon2$WY)
# cre.hydro.wq.mon2=subset(cre.hydro.wq.mon2,WY%in%seq(2010,2018,1))
set.seed(1)
tr.index=sample(1:nrow(cre.hydro.wq.mon),nrow(cre.hydro.wq.mon)*0.7)

cre.hydro.wq.mon[tr.index,"index"]="train"
cre.hydro.wq.mon[-tr.index,"index"]="test"
boxplot(mean.TP~index,cre.hydro.wq.mon)
kruskal.test(mean.TP~index,cre.hydro.wq.mon)

plot(mean.TP~monCY,cre.hydro.wq.mon)
with(cre.hydro.wq.mon[tr.index,],points(monCY,mean.TP,pch=21,bg="red"))

mod=lm(log(mean.TP)~C43Basin+basin.q.ratio+S77+LakeO,cre.hydro.wq.mon[tr.index,])
#mod=lm(log(mean.TP)~grad+C43Basin+S78Qfreq+basin.q.ratio+S79_H,cre.hydro.wq.mon[tr.index,])
layout(matrix(1:4,2,2));plot(mod)
shapiro.test(residuals(mod))
acf(mod$residuals)
vif(mod)
summary(mod)
AIC(mod)

gvlma::gvlma(mod)
lmtest::bgtest(mod)

#broom::glance(mod)
#broom::tidy(mod)
mod.pred=predict(mod,cre.hydro.wq.mon[-tr.index,],interval="confidence")

actuals_preds <-data.frame(cbind(actuals=cre.hydro.wq.mon[-tr.index,"mean.TP"],predicted=exp(mod.pred[,1])))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")


#tiff(filename=paste0(plot.path,"C43TP_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(0.05,0.26);by.x=0.050;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,actuals_preds,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
with(actuals_preds,points(predicted,actuals,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),lwd=0.1,cex=1.25))
mod.cor=mblm(actuals~predicted,actuals_preds)
x.val=with(actuals_preds,seq(min(predicted),max(predicted),length.out = 25))
mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
lines(x.val,mod.cor.pred[,1],col="indianred1",lwd=2)
#lines(x.val,mod.cor.pred[,2],col="indianred1",lty=2)
#lines(x.val,mod.cor.pred[,3],col="indianred1",lty=2)
axis_fun(1,xmaj,xmin,xmaj*1000,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=1,line=1.75,"Predicted TP (\u03BCg L\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TP (\u03BCg L\u207B\u00B9)")
legend("topleft",legend=c("Training Dataset","1:1 line","Correlation"),
       pch=c(21,NA,NA),lty=c(NA,1,1),lwd=c(0.1,1,2),
       col=c(adjustcolor("black",0.5),"black","indianred1"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

#https://rcompanion.org/handbook/G_14.html
#min_max_accuracy
mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  

rcompanion::accuracy(list(mod))
#DescTools::PseudoR2(mod)


dev.off()
plot(mean.TP~monCY,cre.hydro.wq.mon)
with(cre.hydro.wq.mon,pt_line(monCY,mean.TP,2,"grey50",1,21,"dodgerblue1",pt.col="grey50",pt.lwd=0.1))
mod.pred2=predict(mod,cre.hydro.wq.mon,interval="confidence")
shaded.range(cre.hydro.wq.mon$monCY,exp(mod.pred2[,2]),exp(mod.pred2[,3]),"indianred1",lty=1)
lines(cre.hydro.wq.mon$monCY,exp(mod.pred2[,1]),col="indianred1",lwd=1.5)
#lines($monCY,exp(mod.pred2[,2]),col="indianred1",lty=2)
#lines($monCY,exp(mod.pred2[,3]),col="indianred1",lty=2)

#k-folding
#http://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
library(caret)
kmodel=train(
  formula(mod),
  cre.hydro.wq.mon,
  method="lm",
  trControl=trainControl(
    method="cv",
    number=10,
    verboseIter=T
  )
)
kmodel
print(kmodel)
kmodel$results
kmodel$finalModel
summary(kmodel)
shapiro.test(residuals(kmodel))

kfold.mod.pred=predict(kmodel,na.omit(cre.hydro.wq.mon),interval="confidence")
lines(na.omit(cre.hydro.wq.mon)$monCY,exp(kfold.mod.pred),col="black",lwd=2)


kmodel2=data.frame()
kerrors=data.frame()
k=10
set.seed(1)
for(i in 1:k){
  tr.index=sample(1:nrow(cre.hydro.wq.mon),nrow(cre.hydro.wq.mon)*0.7)
  
  train.dat=cre.hydro.wq.mon[tr.index,]
  test.dat=cre.hydro.wq.mon[-tr.index,]
  
  k.cv.mod=lm(log(mean.TP)~C43Basin+basin.q.ratio+S77+LakeO,train.dat)
  summary(k.cv.mod)$r.squared
  
  k.cv.mod.pred=predict(k.cv.mod,test.dat)
  actuals_preds <-data.frame(cbind(actuals=as.numeric(test.dat$mean.TP),predicted=exp(k.cv.mod.pred)))
  MMA.val=mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
  MAPE.val=mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
  
  actuals_preds$k=i
  kmodel2=rbind(kmodel2,actuals_preds)
  kerrors=rbind(kerrors,data.frame(k=i,R2.val=summary(k.cv.mod)$r.squared,RMSE=summary(k.cv.mod)$sigma,MMA=MMA.val,MAPE=MAPE.val))
  print(i)
}
apply(kerrors,2,mean)

ylim.val=c(0,1);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.75,4.25)
plot(1:4,apply(kerrors,2,mean)[2:5],ylim=c(0,1),type="n",axes=F,ylab=NA,xlab=NA,xlim=xlim.val)
arrows(1:4,apply(kerrors[,2:5],2,min),1:4,apply(kerrors[,2:5],2,max),angle=90,length=0.025,code=3,lwd=1)
points(1:4,apply(kerrors,2,mean)[2:5],pch=21,bg="indianred1",lwd=0.1)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:4,1:4,c("R2","RMSE","MMA","MAPE"));box(lwd=1)
mtext(side=2,line=3,"Cross Validation Error")

plot(R2.val~k,kerrors)
plot(MAPE~k,kerrors)
plot(MMA~k,kerrors)

cols=wesanderson::wes_palette("Zissou1",k,"continuous")
#tiff(filename=paste0(plot.path,"C43TP_kmodel.tiff"),width=4,height=2.25,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.25,0.25));
xlim.val=c(0.05,0.26);by.x=0.050;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,kmodel2,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
for(i in 1:k){
  with(subset(kmodel2,k==i),points(predicted,actuals,pch=21,bg=adjustcolor(cols[i],0.5),col=cols[i],lwd=0.1,cex=0.75))
  mod.cor=mblm(actuals~predicted,subset(kmodel2,k==i))
  x.val=with(subset(kmodel2,k==i),seq(min(predicted),max(predicted),length.out = 25))
  mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
  lines(x.val,mod.cor.pred[,1],col=cols[i],lwd=1,lty=2)
}
axis_fun(1,xmaj,xmin,xmaj*1000,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=1,line=1.75,"Predicted TP (\u03BCg L\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TP (\u03BCg L\u207B\u00B9)")
dev.off()


# library(DAAG)
# cv.rslt=CVlm(data=cre.hydro.wq.mon,
#              form.lm= formula(mod),
#              m=10,seed=1,
#              plotit="Observed")
# CVlm(data=cre.hydro.wq.mon,
#              form.lm=formula(mod),
#              m=5,seed=1,
#              plotit="Residual")
# 
# attr(cv.rslt,'ms');# mean square error

# DIN Model ---------------------------------------------------------------
plot(mean.DIN~monCY,cre.hydro.wq.mon)
vars=c("mean.DIN","C43Basin","S79","S77","S78","basin.q.ratio","S78Qfreq","LakeO")
S79.DIN.mod.all=lm(log(mean.DIN)~.,na.omit(cre.hydro.wq.mon[,vars]))
#AIC model
S79.DIN.mod.sw=stepAIC(S79.DIN.mod.all,direction="both",trace=F)
S79.DIN.mod.sw$anova
summary(S79.DIN.mod.sw)
shapiro.test(S79.DIN.mod.sw$residuals)
vif(S79.DIN.mod.sw)

mod.DIN=lm(log(mean.DIN) ~ S78 + basin.q.ratio + S78Qfreq + LakeO ,cre.hydro.wq.mon[tr.index,])
summary(mod.DIN)
layout(matrix(1:4,2,2));plot(mod.DIN)
shapiro.test(residuals(mod.DIN))
vif(mod.DIN)
gvlma::gvlma(mod.DIN)
lmtest::bgtest(mod.DIN)

# TN Model ----------------------------------------------------------------
plot(mean.TN~monCY,cre.hydro.wq.mon)
cre.hydro.wq.mon$S77Q.Qmean=with(cre.hydro.wq.mon,S77/mean(S77,na.rm=T))
cre.hydro.wq.mon$S78Q.Qmean=with(cre.hydro.wq.mon,S78/mean(S78,na.rm=T))

vars=c("mean.TN","C43Basin","S79","S77","S78","basin.q.ratio","S78Qfreq","LakeO")
S79.TN.mod.all=lm(log(mean.TN)~.,na.omit(cre.hydro.wq.mon[,vars]))
layout(matrix(1:4,2,2));plot(S79.TN.mod.all)
shapiro.test(residuals(S79.TN.mod.all))

dev.off()
#AIC model
S79.TN.mod.sw=stepAIC(S79.TN.mod.all,direction="both",trace=F)
S79.TN.mod.sw$anova
summary(S79.TN.mod.sw)
vif(S79.TN.mod.sw)
#StepVIF
library(pedometrics)
S79.TN.mod.sVIF=stepVIF(S79.TN.mod.all)
summary(S79.TN.mod.sVIF)
vif(S79.TN.mod.sVIF)

## Train vs test
boxplot(mean.TN~index,cre.hydro.wq.mon)
kruskal.test(mean.TN~index,cre.hydro.wq.mon)

plot(mean.TN~monCY,cre.hydro.wq.mon)
with(cre.hydro.wq.mon[tr.index,],points(monCY,mean.TN,pch=21,bg="red"))

dev.off()

#mod.TN=lm(log(mean.TN)~S79_H+S77+S78+S78Qfreq,cre.hydro.wq.mon[tr.index,])
mod.TN=lm(log(mean.TN) ~ S79+basin.q.ratio+LakeO ,cre.hydro.wq.mon[tr.index,])
summary(mod.TN)
layout(matrix(1:4,2,2));plot(mod.TN)
shapiro.test(residuals(mod.TN))
acf(mod.TN$residuals)
vif(mod.TN)
AIC(mod.TN)

gvlma::gvlma(mod.TN)
lmtest::bgtest(mod.TN)

mod.TN.pred=predict(mod.TN,cre.hydro.wq.mon[-tr.index,],interval="confidence")

actuals_preds <-data.frame(cbind(actuals=cre.hydro.wq.mon[-tr.index,"mean.TN"],predicted=exp(mod.TN.pred[,1])))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#https://rcompanion.org/handbook/G_14.html
#min_max_accuracy
mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  


#tiff(filename=paste0(plot.path,"C43TN_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(1.1,1.5);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,actuals_preds,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
with(actuals_preds,points(predicted,actuals,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),lwd=0.1,cex=1.25))
mod.cor=mblm(actuals~predicted,actuals_preds)
x.val=with(actuals_preds,seq(min(predicted),max(predicted),length.out = 25))
mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
lines(x.val,mod.cor.pred[,1],col="indianred1",lwd=2)
#lines(x.val,mod.cor.pred[,2],col="indianred1",lty=2)
#lines(x.val,mod.cor.pred[,3],col="indianred1",lty=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=1.75,"Predicted TN (mg L\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TN (mg L\u207B\u00B9)")
legend("topleft",legend=c("Training Dataset","1:1 line","Correlation"),
       pch=c(21,NA,NA),lty=c(NA,1,1),lwd=c(0.1,1,2),
       col=c(adjustcolor("black",0.5),"black","indianred1"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

dev.off()
plot(mean.TN~monCY,cre.hydro.wq.mon)
with(cre.hydro.wq.mon,pt_line(monCY,mean.TN,2,"grey50",1,21,"dodgerblue1",pt.col="grey50",pt.lwd=0.1))
mod.pred2=predict(mod.TN,cre.hydro.wq.mon,interval="confidence")
shaded.range(cre.hydro.wq.mon$monCY,exp(mod.pred2[,2]),exp(mod.pred2[,3]),"indianred1",lty=1)
lines(cre.hydro.wq.mon$monCY,exp(mod.pred2[,1]),col="indianred1",lwd=1.5)

library(caret)
kmodel=train(
  formula(S79.TN.mod.sw),
  na.omit(cre.hydro.wq.mon),
  method="lm",
  trControl=trainControl(
    method="cv",
    number=10,
    verboseIter=T
  )
)
kmodel
print(kmodel)
kmodel$results

plot(mean.TN~monCY,cre.hydro.wq.mon)
with(cre.hydro.wq.mon,pt_line(monCY,mean.TN,2,"grey50",1,21,"dodgerblue1",pt.col="grey50",pt.lwd=0.1))
mod.pred2=predict(mod.TN,cre.hydro.wq.mon,interval="confidence")
shaded.range(cre.hydro.wq.mon$monCY,exp(mod.pred2[,2]),exp(mod.pred2[,3]),"indianred1",lty=1)
lines(cre.hydro.wq.mon$monCY,exp(mod.pred2[,1]),col="indianred1",lwd=1.5)
kfold.mod.pred=predict(kmodel,na.omit(cre.hydro.wq.mon),interval="confidence")
lines(na.omit(cre.hydro.wq.mon)$monCY,exp(kfold.mod.pred),col="black",lwd=1.5)


kmodel2.TN=data.frame()
TN.kerrors=data.frame()
k=10
set.seed(1)
for(i in 1:k){
  tr.index=sample(1:nrow(cre.hydro.wq.mon),nrow(cre.hydro.wq.mon)*0.7)
  
  train.dat=cre.hydro.wq.mon[tr.index,]
  test.dat=cre.hydro.wq.mon[-tr.index,]
  
  k.cv.mod=lm(log(mean.TN) ~ S79+basin.q.ratio+LakeO,train.dat)
  
  k.cv.mod.pred=predict(k.cv.mod,test.dat)
  actuals_preds <-data.frame(cbind(actuals=as.numeric(test.dat$mean.TN),predicted=exp(k.cv.mod.pred)))
  MMA.val=mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
  MAPE.val=mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
  
  actuals_preds$k=i
  kmodel2.TN=rbind(kmodel2.TN,actuals_preds)
  TN.kerrors=rbind(TN.kerrors,data.frame(k=i,R2.val=summary(k.cv.mod)$r.squared,RMSE=summary(k.cv.mod)$sigma,MMA=MMA.val,MAPE=MAPE.val))
  print(i)
}
apply(TN.kerrors,2,mean)

ylim.val=c(0,1);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.75,4.25)
plot(1:4,apply(TN.kerrors,2,mean)[2:5],ylim=c(0,1),type="n",axes=F,ylab=NA,xlab=NA,xlim=xlim.val)
arrows(1:4,apply(TN.kerrors[,2:5],2,min),1:4,apply(TN.kerrors[,2:5],2,max),angle=90,length=0.025,code=3,lwd=1)
points(1:4,apply(TN.kerrors,2,mean)[2:5],pch=21,bg="indianred1",lwd=0.1)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:4,1:4,c("R2","RMSE","MMA","MAPE"));box(lwd=1)
mtext(side=2,line=3,"Cross Validation Error")

plot(R2.val~k,TN.kerrors)
plot(MAPE~k,TN.kerrors)
plot(MMA~k,TN.kerrors)
plot(actuals~predicted,kmodel2.TN)

cols=wesanderson::wes_palette("Zissou1",k,"continuous")
#tiff(filename=paste0(plot.path,"C43TN_kmodel.tiff"),width=4,height=2.25,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.25,0.25));
xlim.val=c(1.1,1.6);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0.9,1.6);by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,kmodel2.TN,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
for(i in 1:k){
  with(subset(kmodel2.TN,k==i),points(predicted,actuals,pch=21,bg=adjustcolor(cols[i],0.5),col=cols[i],lwd=0.1,cex=0.75))
  mod.cor=mblm(actuals~predicted,subset(kmodel2.TN,k==i))
  x.val=with(subset(kmodel2.TN,k==i),seq(min(predicted),max(predicted),length.out = 25))
  mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
  lines(x.val,mod.cor.pred[,1],col=cols[i],lwd=1,lty=2)
}
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.75,"Predicted TN (mg L\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TN (mg L\u207B\u00B9)")
dev.off()


# TON Model ---------------------------------------------------------------
plot(mean.TON~monCY,cre.hydro.wq.mon)
vars=c("mean.TON","C43Basin","S79","S77","S78","basin.q.ratio","S78Qfreq","LakeO")
S79.TON.mod.all=lm(log(mean.TON)~.,na.omit(cre.hydro.wq.mon[,vars]))
#AIC model
S79.TON.mod.sw=stepAIC(S79.TON.mod.all,direction="both",trace=F)
S79.TON.mod.sw$anova
summary(S79.TON.mod.sw)
shapiro.test(S79.TON.mod.sw$residuals)
vif(S79.TON.mod.sw)

mod.TON=lm(log(mean.TON) ~ S78 + C43Basin+LakeO ,cre.hydro.wq.mon[tr.index,])
summary(mod.TON)
layout(matrix(1:4,2,2));plot(mod.TON)
shapiro.test(residuals(mod.TON))
vif(mod.TON)
gvlma::gvlma(mod.TON)
lmtest::bgtest(mod.TON)

mod.TON.pred=predict(mod.TON,cre.hydro.wq.mon[-tr.index,],interval="confidence")
actuals_preds <-data.frame(cbind(actuals=cre.hydro.wq.mon[-tr.index,"mean.TON"],predicted=exp(mod.TON.pred[,1])))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")
#https://rcompanion.org/handbook/G_14.html
#min_max_accuracy
mean(apply(actuals_preds, 1, min,na.rm=T) / apply(actuals_preds, 1, max,na.rm=T)) 
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals,na.rm=T)  

#tiff(filename=paste0(plot.path,"C43TON_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(0.8,1.5);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,actuals_preds,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
with(actuals_preds,points(predicted,actuals,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),lwd=0.1,cex=1.25))
mod.cor=mblm(actuals~predicted,na.omit(actuals_preds))
x.val=with(actuals_preds,seq(min(predicted),max(predicted),length.out = 25))
mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
lines(x.val,mod.cor.pred[,1],col="indianred1",lwd=2)
#lines(x.val,mod.cor.pred[,2],col="indianred1",lty=2)
#lines(x.val,mod.cor.pred[,3],col="indianred1",lty=2)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.75,"Predicted TON (mg L\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TON (mg L\u207B\u00B9)")
legend("topleft",legend=c("Training Dataset","1:1 line","Correlation"),
       pch=c(21,NA,NA),lty=c(NA,1,1),lwd=c(0.1,1,2),
       col=c(adjustcolor("black",0.5),"black","indianred1"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

kmodel2.TON=data.frame()
TON.kerrors=data.frame()
k=10
set.seed(1)
for(i in 1:k){
  tr.index=sample(1:nrow(cre.hydro.wq.mon),nrow(cre.hydro.wq.mon)*0.7)
  
  train.dat=cre.hydro.wq.mon[tr.index,]
  test.dat=cre.hydro.wq.mon[-tr.index,]
  
  k.cv.mod=lm(log(mean.TON) ~ S78 + C43Basin+LakeO ,train.dat)
  
  k.cv.mod.pred=predict(k.cv.mod,test.dat)
  actuals_preds <-na.omit(data.frame(cbind(actuals=as.numeric(test.dat$mean.TON),predicted=exp(k.cv.mod.pred))))
  MMA.val=mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
  MAPE.val=mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
  
  actuals_preds$k=i
  kmodel2.TON=rbind(kmodel2.TON,actuals_preds)
  TON.kerrors=rbind(TON.kerrors,data.frame(k=i,R2.val=summary(k.cv.mod)$r.squared,RMSE=summary(k.cv.mod)$sigma,MMA=MMA.val,MAPE=MAPE.val))
  print(i)
}
apply(TON.kerrors,2,mean)

ylim.val=c(0,1);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.75,4.25)
plot(1:4,apply(TON.kerrors,2,mean)[2:5],ylim=c(0,1),type="n",axes=F,ylab=NA,xlab=NA,xlim=xlim.val)
arrows(1:4,apply(TON.kerrors[,2:5],2,min),1:4,apply(TON.kerrors[,2:5],2,max),angle=90,length=0.025,code=3,lwd=1)
points(1:4,apply(TON.kerrors,2,mean)[2:5],pch=21,bg="indianred1",lwd=0.1)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:4,1:4,c("R2","RMSE","MMA","MAPE"));box(lwd=1)
mtext(side=2,line=3,"Cross Validation Error")

plot(R2.val~k,TON.kerrors)
plot(MAPE~k,TON.kerrors)
plot(MMA~k,TON.kerrors)
plot(actuals~predicted,kmodel2.TON)

cols=wesanderson::wes_palette("Zissou1",k,"continuous")
#tiff(filename=paste0(plot.path,"C43TON_kmodel.tiff"),width=4,height=2.25,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.25,0.25));
xlim.val=c(0.8,1.5);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,kmodel2.TON,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
for(i in 1:k){
  with(subset(kmodel2.TON,k==i),points(predicted,actuals,pch=21,bg=adjustcolor(cols[i],0.5),col=cols[i],lwd=0.1,cex=0.75))
  mod.cor=mblm(actuals~predicted,na.omit(subset(kmodel2.TON,k==i)))
  x.val=with(subset(kmodel2.TON,k==i),seq(min(predicted),max(predicted),length.out = 25))
  mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
  lines(x.val,mod.cor.pred[,1],col=cols[i],lwd=1,lty=2)
}
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.75,"Predicted TON (mg L\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TON (mg L\u207B\u00B9)")
dev.off()


# Load Est ----------------------------------------------------------------
wq.dat.xtab.mon2=ddply(subset(wq.dat.xtab,Station.ID=="S79"),c("WY","monCY"),summarise,mean.TP=mean(TP,na.rm=T),mean.TN=mean(TN,na.rm=T),mean.DIN=mean(DIN,na.rm=T),mean.TON=mean(TON,na.rm=T),mean.TNTP=mean(TNTP,na.rm=T))
cre.hydro.wq.mon2=merge(cre.hydro.mon,wq.dat.xtab.mon2,c("WY","monCY"),all.x=T)
diff(cre.hydro.wq.mon2$monCY)
range(cre.hydro.wq.mon2$WY)
cre.hydro.wq.mon2=subset(cre.hydro.wq.mon2,WY>2009)

cre.hydro.wq.mon2$TP.predict=exp(predict(mod,cre.hydro.wq.mon2,interval="confidence"))
cre.hydro.wq.mon2$Obs.TPLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,mean.TP))
cre.hydro.wq.mon2$Pred.fit.TPLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TP.predict)[,1]))
cre.hydro.wq.mon2$Pred.LCI.TPLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TP.predict)[,2]))
cre.hydro.wq.mon2$Pred.UCI.TPLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TP.predict)[,3]))

cre.hydro.wq.mon2$TN.predict=exp(predict(mod.TN,cre.hydro.wq.mon2,interval="confidence"))
cre.hydro.wq.mon2$Obs.TNLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,mean.TN))
cre.hydro.wq.mon2$Pred.fit.TNLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TN.predict)[,1]))
cre.hydro.wq.mon2$Pred.LCI.TNLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TN.predict)[,2]))
cre.hydro.wq.mon2$Pred.UCI.TNLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TN.predict)[,3]))

cre.hydro.wq.mon2$TON.predict=exp(predict(mod.TON,cre.hydro.wq.mon2,interval="confidence"))
cre.hydro.wq.mon2$Obs.TONLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,mean.TON))
cre.hydro.wq.mon2$Pred.fit.TONLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TON.predict)[,1]))
cre.hydro.wq.mon2$Pred.LCI.TONLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TON.predict)[,2]))
cre.hydro.wq.mon2$Pred.UCI.TONLoad.kg=with(cre.hydro.wq.mon2,Load.Calc.kg(S79/2.44658e-6,unlist(TON.predict)[,3]))

with(cre.hydro.wq.mon2,cor.test(Pred.fit.TPLoad.kg,Obs.TPLoad.kg,method="spearman"))
with(cre.hydro.wq.mon2,cor.test(Pred.fit.TNLoad.kg,Obs.TNLoad.kg,method="spearman"))
with(cre.hydro.wq.mon2,cor.test(Pred.fit.TONLoad.kg,Obs.TONLoad.kg,method="spearman"))

#tiff(filename=paste0(plot.path,"C43_ObsPredloads.tiff"),width=6.5,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/C43_ObsPredloads.png"),width=6.5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.25,1),oma=c(2,2.5,0.75,0.25));
layout(matrix(1:3,1,3))
txt.cex=0.8
xlim.val=c(0,46)*1000;by.x=10000;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.fit.TPLoad.kg~Obs.TPLoad.kg,cre.hydro.wq.mon2,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1.25,lty=2)
with(cre.hydro.wq.mon2,arrows(Obs.TPLoad.kg,Pred.LCI.TPLoad.kg,Obs.TPLoad.kg,Pred.UCI.TPLoad.kg,angle=90,code=3,length=0,col=adjustcolor("black",0.25)))
with(cre.hydro.wq.mon2,points(Obs.TPLoad.kg,Pred.fit.TPLoad.kg,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),cex=1.25,lwd=0.1))
axis_fun(1,xmaj,xmin,format(xmaj/1000),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1000));box(lwd=1)
mtext(side=2,line=2,"Pred. TP Load (x10\u00B3 kg mon\u207B\u00B9)",cex=txt.cex)
mtext(side=1,line=1.75,"Obs. TP Load (x10\u00B3 kg mon\u207B\u00B9)",cex=txt.cex)

xlim.val=c(0,12)*100000;by.x=300000;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.fit.TNLoad.kg~Obs.TNLoad.kg,cre.hydro.wq.mon2,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1.25,lty=2)
with(cre.hydro.wq.mon2,arrows(Obs.TNLoad.kg,Pred.LCI.TNLoad.kg,Obs.TNLoad.kg,Pred.UCI.TNLoad.kg,angle=90,code=3,length=0,col=adjustcolor("black",0.25)))
with(cre.hydro.wq.mon2,points(Obs.TNLoad.kg,Pred.fit.TNLoad.kg,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),cex=1.25,lwd=0.1))
axis_fun(1,xmaj,xmin,format(xmaj/1e5),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=2,line=2,"Pred. TN Load (x10\u2074 kg mon\u207B\u00B9)",cex=txt.cex)
mtext(side=1,line=1.75,"Obs. TN Load (x10\u2074 kg mon\u207B\u00B9)",cex=txt.cex)

plot(Pred.fit.TONLoad.kg~Obs.TONLoad.kg,cre.hydro.wq.mon2,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1.25,lty=2)
with(cre.hydro.wq.mon2,arrows(Obs.TONLoad.kg,Pred.LCI.TONLoad.kg,Obs.TONLoad.kg,Pred.UCI.TONLoad.kg,angle=90,code=3,length=0,col=adjustcolor("black",0.25)))
with(cre.hydro.wq.mon2,points(Obs.TONLoad.kg,Pred.fit.TONLoad.kg,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),cex=1.25,lwd=0.1))
axis_fun(1,xmaj,xmin,format(xmaj/1e5),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=2,line=2,"Pred. TON Load (x10\u2074 kg mon\u207B\u00B9)",cex=txt.cex)
mtext(side=1,line=1.75,"Obs. TON Load (x10\u2074 kg mon\u207B\u00B9)",cex=txt.cex)
dev.off()


WYLoads=ddply(cre.hydro.wq.mon2,"WY",summarise,
              Obs.TP=sum(Obs.TPLoad.kg,na.rm=T),Pred.TP=sum(Pred.fit.TPLoad.kg,na.rm=T),LCI.TP=sum(Pred.LCI.TPLoad.kg,na.rm=T),UCI.TP=sum(Pred.UCI.TPLoad.kg,na.rm=T),
              Obs.TN=sum(Obs.TNLoad.kg,na.rm=T),Pred.TN=sum(Pred.fit.TNLoad.kg,na.rm=T),LCI.TN=sum(Pred.LCI.TNLoad.kg,na.rm=T),UCI.TN=sum(Pred.UCI.TNLoad.kg,na.rm=T),
              Obs.TON=sum(Obs.TONLoad.kg,na.rm=T),Pred.TON=sum(Pred.fit.TONLoad.kg,na.rm=T),LCI.TON=sum(Pred.LCI.TONLoad.kg,na.rm=T),UCI.TON=sum(Pred.UCI.TONLoad.kg,na.rm=T))

summary(WYLoads)

#tiff(filename=paste0(plot.path,"C43_WYObsPredloads.tiff"),width=6.5,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,0.25,1),oma=c(2,1.5,0.75,0.25));
layout(matrix(1:3,1,3))
txt.cex=0.8

xlim.val=c(2010,2019);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,60)*10000;by.y=200000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Obs.TP~WY,WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(WYLoads,shaded.range(WY,LCI.TP,UCI.TP,"forestgreen",lty=1))
with(WYLoads,pt_line(WY,Pred.TP,2,"forestgreen",1,21,"forestgreen"))
with(WYLoads,pt_line(WY,Obs.TP,2,"indianred1",1,21,"indianred1"))
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e4));box(lwd=1)
mtext(side=2,line=2,"TP Load (x10\u00B3 kg WY\u207B\u00B9)",cex=txt.cex)
#mtext(side=1,line=1.75,"Florida Water Year (May - Apirl)",cex=txt.cex)
legend("topleft",legend=c("Observed Load","Predicted Load \u00B1 95% CI"),
       pch=21,lty=NA,lwd=0.1,
       col="black",pt.bg=c("indianred1","forestgreen"),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(0,50)*100000;by.y=1000000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Obs.TN~WY,WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(WYLoads,shaded.range(WY,LCI.TN,UCI.TN,"forestgreen",lty=1))
with(WYLoads,pt_line(WY,Pred.TN,2,"forestgreen",1,21,"forestgreen"))
with(WYLoads,pt_line(WY,Obs.TN,2,"indianred1",1,21,"indianred1"))
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=2,line=2,"TN Load (x10\u2074 kg WY\u207B\u00B9)",cex=txt.cex)
mtext(side=1,line=1.75,"Florida Water Year (May - Apirl)",cex=txt.cex)

plot(Obs.TON~WY,WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(WYLoads,shaded.range(WY,LCI.TON,UCI.TON,"forestgreen",lty=1))
with(WYLoads,pt_line(WY,Pred.TON,2,"forestgreen",1,21,"forestgreen"))
with(WYLoads,pt_line(WY,Obs.TON,2,"indianred1",1,21,"indianred1"))
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=2,line=2,"TON Load (x10\u2074 kg WY\u207B\u00B9)",cex=txt.cex)
#mtext(side=1,line=1.75,"Florida Water Year (May - Apirl)",cex=txt.cex)
dev.off()


# RSM  --------------------------------------------------------------------
library(dssrip)

dss_out=opendss(paste0(data.path,"RSMBN/pro_19652016_LORS08/RSMBN_output.dss"))
paths=paste0("/RSMBN/LOK/STAGE/01JAN1965 - 01JAN2016/1DAY/SIMULATED/")
RSM.lakeO.stg=data.frame(getFullTSC(dss_out,paths))
RSM.lakeO.stg$Date=date.fun(rownames(RSM.lakeO.stg))

RSM.lakeO.stg$WY=WY(RSM.lakeO.stg$Date)
RSM.lakeO.stg$month=format(RSM.lakeO.stg$Date,"%m")
RSM.lakeO.stg$CY=format(RSM.lakeO.stg$Date,"%Y")
RSM.lakeO.stg$monCY=with(RSM.lakeO.stg,date.fun(paste(CY,month,"01",sep="-")))

RSM.lakeO.stg.mon=ddply(RSM.lakeO.stg,c("monCY","WY"),summarise,
                    LakeO=mean(STAGE,na.rm=T))
plot(LakeO~monCY,RSM.lakeO.stg.mon)

RSM.sites=c("S77","S78","S79")
RSM.q.cre.dat=data.frame()
for(i in 1:length(RSM.sites)){
  paths=paste0("/RSMBN/",RSM.sites[i],"/FLOW/01JAN1965 - 01JAN2016/1DAY/SIMULATED/")  
  tmp=data.frame(getFullTSC(dss_out,paths))
  tmp$Date=date.fun(rownames(RSM.lakeO.stg))
  tmp$SITE=RSM.sites[i]
  RSM.q.cre.dat=rbind(tmp,RSM.q.cre.dat)
  print(i)
}

RSM.q.cre.dat.xtab=cast(RSM.q.cre.dat,Date~SITE,value="FLOW",mean)
RSM.q.cre.dat.xtab$month=format(RSM.q.cre.dat.xtab$Date,"%m")
RSM.q.cre.dat.xtab$CY=format(RSM.q.cre.dat.xtab$Date,"%Y")
RSM.q.cre.dat.xtab$monCY=with(RSM.q.cre.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
RSM.q.cre.dat.xtab$WY=WY(RSM.q.cre.dat.xtab$Date)
RSM.q.cre.dat.xtab$S78QFreq=as.numeric(RSM.q.cre.dat.xtab$S78>0)
summary(RSM.q.cre.dat.xtab)

RSM.q.cre.dat.xtab.mon=ddply(RSM.q.cre.dat.xtab,c("monCY","WY"),summarise,S77=sum(cfs.to.km3d(S77),na.rm=T),S78=sum(cfs.to.km3d(S78),na.rm=T),S79=sum(cfs.to.km3d(S79),na.rm=T),S78Qfreq=sum(S78QFreq,na.rm=T))
RSM.q.cre.dat.xtab.mon$C43Basin=with(RSM.q.cre.dat.xtab.mon,ifelse(S79<S77,0,S79-S77))
RSM.q.cre.dat.xtab.mon$LakeQ=with(RSM.q.cre.dat.xtab.mon,ifelse(S79-S77>0,S77,S79))
RSM.q.cre.dat.xtab.mon$basin.q.ratio=with(RSM.q.cre.dat.xtab.mon,ifelse(S79==0,0,C43Basin/S79))

plot(S79~monCY,RSM.q.cre.dat.xtab.mon)
plot(C43Basin~monCY,RSM.q.cre.dat.xtab.mon)

RSM.cre.hydro=merge(RSM.q.cre.dat.xtab.mon,RSM.lakeO.stg.mon,c("monCY","WY"))
ddply(RSM.cre.hydro,"WY",summarise,N.val=N.obs(S79))
RSM.cre.hydro=subset(RSM.cre.hydro,WY%in%seq(1966,2016,1))

Load.Calc.kg.LOSOM=function(flow.km3d,conc.mgL){
  q.cmd=flow.km3d*1e9
  conc.mgm3=conc.mgL*1000
  load.kg=(q.cmd*conc.mgm3)*1e-6
  return(load.kg)
}

RSM.cre.hydro$TP.fit=exp(predict(mod,RSM.cre.hydro,interval="confidence"))[,1]
RSM.cre.hydro$TP.LCI=exp(predict(mod,RSM.cre.hydro,interval="confidence"))[,2]
RSM.cre.hydro$TP.UCI=exp(predict(mod,RSM.cre.hydro,interval="confidence"))[,3]
RSM.cre.hydro$Pred.fit.TPLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TP.fit))
RSM.cre.hydro$Pred.LCI.TPLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TP.LCI))
RSM.cre.hydro$Pred.UCI.TPLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TP.UCI))

RSM.cre.hydro$TN.fit=exp(predict(mod.TN,RSM.cre.hydro,interval="confidence"))[,1]
RSM.cre.hydro$TN.LCI=exp(predict(mod.TN,RSM.cre.hydro,interval="confidence"))[,2]
RSM.cre.hydro$TN.UCI=exp(predict(mod.TN,RSM.cre.hydro,interval="confidence"))[,3]
RSM.cre.hydro$Pred.fit.TNLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TN.fit))
RSM.cre.hydro$Pred.LCI.TNLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TN.LCI))
RSM.cre.hydro$Pred.UCI.TNLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TN.UCI))

RSM.cre.hydro$TON.fit=exp(predict(mod.TON,RSM.cre.hydro,interval="confidence"))[,1]
RSM.cre.hydro$TON.LCI=exp(predict(mod.TON,RSM.cre.hydro,interval="confidence"))[,2]
RSM.cre.hydro$TON.UCI=exp(predict(mod.TON,RSM.cre.hydro,interval="confidence"))[,3]
RSM.cre.hydro$Pred.fit.TONLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TON.fit))
RSM.cre.hydro$Pred.LCI.TONLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TON.LCI))
RSM.cre.hydro$Pred.UCI.TONLoad.kg=with(RSM.cre.hydro,Load.Calc.kg.LOSOM(S79,TON.UCI))


RSM.WYLoads=ddply(RSM.cre.hydro,"WY",summarise,
              TFlow.km3Y=sum(S79,na.rm=T),
              Pred.TP=sum(Pred.fit.TPLoad.kg,na.rm=T),LCI.TP=sum(Pred.LCI.TPLoad.kg,na.rm=T),UCI.TP=sum(Pred.UCI.TPLoad.kg,na.rm=T),
              Pred.TN=sum(Pred.fit.TNLoad.kg,na.rm=T),LCI.TN=sum(Pred.LCI.TNLoad.kg,na.rm=T),UCI.TN=sum(Pred.UCI.TNLoad.kg,na.rm=T),
              Pred.TON=sum(Pred.fit.TONLoad.kg,na.rm=T),LCI.TON=sum(Pred.LCI.TONLoad.kg,na.rm=T),UCI.TON=sum(Pred.UCI.TONLoad.kg,na.rm=T))

RSM.WYLoads$Pred.TPFWM.ugL=with(RSM.WYLoads,(Pred.TP*1000000)/(TFlow.km3Y*1e9))
RSM.WYLoads$Pred.TNFWM.mgL=with(RSM.WYLoads,(Pred.TN*1000000)/(TFlow.km3Y*1e9))*0.001
RSM.WYLoads$Pred.TONFWM.mgL=with(RSM.WYLoads,(Pred.TON*1000000)/(TFlow.km3Y*1e9))*0.001

#tiff(filename=paste0(plot.path,"RSM_C43_WYPredloads.tiff"),width=6.5,height=5.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,0.25,1),oma=c(2,2,0.75,0.25));
layout(matrix(1:6,3,2,byrow = T))
txt.cex=0.8

xlim.val=c(1966,2016);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,65)*10000;by.y=200000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.TP~WY,RSM.WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.WYLoads,pt_line(WY,Pred.TP,2,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e4));box(lwd=1)
mtext(side=2,line=2.5,"TP Load\n(x10\u00B3 kg WY\u207B\u00B9)",cex=txt.cex)
legend("topleft",legend=c("LOSOM BN - LORS08"),
       pch=21,lty=NA,lwd=0.1,
       col="black",pt.bg=c("forestgreen"),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(100,200);by.y=25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.TPFWM.ugL~WY,RSM.WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.WYLoads,pt_line(WY,Pred.TPFWM.ugL,2,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"TP FWM (\u03BCg L\u207B\u00B9)",cex=txt.cex)

ylim.val=c(0,60)*100000;by.y=2000000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.TN~WY,RSM.WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.WYLoads,pt_line(WY,Pred.TN,2,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e4));box(lwd=1)
mtext(side=2,line=2.5,"TN Load\n(x10\u2074 kg WY\u207B\u00B9)",cex=txt.cex)

ylim.val=c(1,1.4);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.TNFWM.mgL~WY,RSM.WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.WYLoads,pt_line(WY,Pred.TNFWM.mgL,2,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"TN FWM (mg L\u207B\u00B9)",cex=txt.cex)

ylim.val=c(0,60)*100000;by.y=2000000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.TON~WY,RSM.WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.WYLoads,pt_line(WY,Pred.TON,2,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e4));box(lwd=1)
mtext(side=2,line=2.5,"TON Load\n(x10\u2074 kg WY\u207B\u00B9)",cex=txt.cex)

ylim.val=c(1,1.4);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.TONFWM.mgL~WY,RSM.WYLoads,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.WYLoads,pt_line(WY,Pred.TONFWM.mgL,2,"forestgreen",1,21,"forestgreen"))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"TON FWM (mg L\u207B\u00B9)",cex=txt.cex)
mtext(side=1,"Florida Water Year",outer=T,line=1)
dev.off()

# S79 hind cast models ----------------------------------------------------
### hind cast eval
# hind.hydro=subset(cre.hydro.mon,WY<2010)
# hind.wq=subset(wq.dat.xtab.mon,WY<2010)
# hind.wq.hydro=merge(hind.wq,hind.hydro,c("WY","monCY"))
# 
# hind.TN=predict(mod.TN,hind.wq.hydro,interval="confidence")
# plot(mean.TN~monCY,hind.wq.hydro)
# points(hind.wq.hydro$monCY,exp(hind.TN[,1]),pch=21,bg="dodgerblue1")
# mean(apply(cbind(hind.wq$mean.TN,exp(hind.TN[,1])), 1, min,na.rm=T) / apply(cbind(hind.wq$mean.TN,exp(hind.TN[,1])), 1, max,na.rm=T)) 
# mean(abs((exp(hind.TN[,1]) - hind.wq$mean.TN))/hind.wq$mean.TN,na.rm=T)  
# 
# hind.TP=predict(mod,hind.wq.hydro,interval="confidence")
# plot(mean.TP~monCY,hind.wq.hydro,type="b")
# shaded.range(hind.wq.hydro$monCY,exp(hind.TP[,2]),exp(hind.TP[,3]),"grey")
# points(hind.wq.hydro$monCY,exp(hind.TP[,1]),pch=21,bg="dodgerblue1")
# mean(apply(cbind(hind.wq$mean.TP,exp(hind.TP[,1])), 1, min,na.rm=T) / apply(cbind(hind.wq$mean.TP,exp(hind.TP[,1])), 1, max,na.rm=T)) 
# mean(abs((exp(hind.TP[,1]) - hind.wq$mean.TP))/hind.wq$mean.TP,na.rm=T)  
# 
