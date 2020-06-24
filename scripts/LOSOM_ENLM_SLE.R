## 
## LOSOM Water Quality Subteam
## Estuary Nutrient Load Model
## St Lucie Estuary (S80)
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
# SLE ---------------------------------------------------------------------
# Discharge ---------------------------------------------------------------
q.dbkeys=data.frame(SITE=c("S308","S80"),DBKEY=c("15626","JW224"))
q.sle.dat=data.frame()
for(i in 1:nrow(q.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],q.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(q.dbkeys$DBKEY[i])
  q.sle.dat=rbind(q.sle.dat,tmp)
  print(i)
}
q.sle.dat=merge(q.sle.dat,q.dbkeys,"DBKEY")
q.sle.dat$Date.EST=date.fun(q.sle.dat$Date)

#range(q.sle.dat$Data.Value,na.rm=T)
q.sle.dat$Data.Value=with(q.sle.dat,ifelse(Data.Value<0,0,as.numeric(Data.Value)))

q.sle.dat.xtab=cast(q.sle.dat,Date.EST~SITE,value="Data.Value",mean)
q.sle.dat.xtab$month=format(q.sle.dat.xtab$Date,"%m")
q.sle.dat.xtab$CY=format(q.sle.dat.xtab$Date,"%Y")
q.sle.dat.xtab$monCY=with(q.sle.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
q.sle.dat.xtab$WY=WY(q.sle.dat.xtab$Date.EST)
q.sle.dat.xtab$C44.basin.in=with(q.sle.dat.xtab,ifelse(S80<S308,0,S80-S308))
q.sle.dat.xtab$basin.ratio=with(q.sle.dat.xtab,ifelse(S80==0,NA,C44.basin.in/S80))

plot(S308~Date.EST,q.sle.dat.xtab)
plot(S80~Date.EST,q.sle.dat.xtab)
plot(C44.basin.in~Date.EST,q.sle.dat.xtab)
plot(basin.ratio~Date.EST,q.sle.dat.xtab)

# Monthly aggregation
q.sle.dat.xtab.mon=ddply(q.sle.dat.xtab,c("monCY","WY"),summarise,S308=sum(cfs.to.km3d(S308),na.rm=T),S80=sum(cfs.to.km3d(S80),na.rm=T))
q.sle.dat.xtab.mon$C44Basin=with(q.sle.dat.xtab.mon,ifelse(S80<S308,0,S80-S308))
q.sle.dat.xtab.mon$basin.q.ratio=with(q.sle.dat.xtab.mon,C44Basin/S80)


# Lake O Stage ------------------------------------------------------------
lakeo.stg.dbkey=data.frame(SITE=c("LakeO","LakeO"),DBKEY=c("15611","06832"),priority=c("P1","P2"));#15611 Pref DBKEY
lakeO.stg=DBHYDRO_daily(dates[1],dates[2],lakeo.stg.dbkey$DBKEY[1])
plot(Data.Value~Date,subset(lakeO.stg,DBKEY=="15611"))

lakeO.stg$Date.EST=date.fun(lakeO.stg$Date)

lakeO.stg$WY=WY(lakeO.stg$Date)
lakeO.stg$month=format(lakeO.stg$Date,"%m")
lakeO.stg$CY=format(lakeO.stg$Date,"%Y")
lakeO.stg$monCY=with(lakeO.stg,date.fun(paste(CY,month,"01",sep="-")))
lakeO.stg$delta.stg=with(lakeO.stg,c(NA,diff(Data.Value)))
# Monthly aggregation
lakeO.stg.mon=ddply(lakeO.stg,c("monCY","WY"),summarise,
                    LakeO=mean(Data.Value,na.rm=T),
                    del.stg=mean(delta.stg,na.rm=T))
plot(LakeO~monCY,lakeO.stg.mon)

# Water Quality -----------------------------------------------------------
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7),param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp"))
wq.sites=c("C44S80","S308C")
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
#tiff(filename=paste0(plot.path,"S80_wq_monitoring.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/S80_wq_monitoring.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
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
unique(q.sle.dat$SITE)
q.sle.dat$Q.cfs=q.sle.dat$Data.Value

wq.dat.xtab=merge(wq.dat.xtab,data.frame(Station.ID=wq.sites,SITE=c("S80","S308")),"Station.ID")
wq.dat.xtab=merge(wq.dat.xtab,q.sle.dat[,c("SITE","Date.EST","Q.cfs")],c("SITE","Date.EST"),all.x=T)
range(wq.dat.xtab$Q.cfs,na.rm=T)

# Reversal Evaluation
wq.dat.xtab$TPReversal=with(wq.dat.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
wq.dat.xtab$TNReversal=with(wq.dat.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

sum(wq.dat.xtab$TNReversal,na.rm=T)
sum(wq.dat.xtab$TPReversal,na.rm=T)

range(wq.dat.xtab$TN,na.rm=T)
range(wq.dat.xtab$TP,na.rm=T)
subset(wq.dat.xtab,TN>5)
subset(wq.dat.xtab,TP>0.8)

wq.dat.xtab$TN=with(wq.dat.xtab,ifelse(TN>6,NA,TN))

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,wq.dat.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,wq.dat.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

dev.off()

# Monthly aggregation
wq.dat.xtab.mon=ddply(subset(wq.dat.xtab,Station.ID=="C44S80"&abs(Q.cfs)>0),c("WY","monCY"),summarise,mean.TP=mean(TP,na.rm=T),mean.TN=mean(TN,na.rm=T),mean.DIN=mean(DIN,na.rm=T),mean.TON=mean(TON,na.rm=T))
#wq.dat.xtab.mon=ddply(subset(wq.dat.xtab,Station.ID=="C44S80"),c("WY","monCY"),summarise,mean.TP=mean(TP,na.rm=T),mean.TN=mean(TN,na.rm=T),mean.DIN=mean(DIN,na.rm=T),mean.TON=mean(TON,na.rm=T))
wq.dat.xtab.mon

# Concentration Discharge Relationship ------------------------------------
mean(subset(q.sle.dat.xtab,WY>2009)$S80,na.rm=T)
mean(q.sle.dat.xtab$S80,na.rm=T)
q.sle.dat.xtab$S80.QQmean=with(q.sle.dat.xtab,S80/mean(subset(q.sle.dat.xtab,WY>2009)$S80,na.rm=T))
q.sle.dat.xtab$C44.QQmean=with(q.sle.dat.xtab,C44.basin.in/mean(subset(q.sle.dat.xtab,WY>2009)$C44.basin.in,na.rm=T))

plot(S80~Date.EST,q.sle.dat.xtab,type="l")
abline(h=mean(subset(q.sle.dat.xtab,WY>2009)$S80,na.rm=T),lty=2,col="indianred1")
abline(h=mean(q.sle.dat.xtab$S80,na.rm=T),lty=2,col="indianred1")

wq.vars=c("Station.ID","WY", "Date.EST", "Chla", "NH4", "NOx", "SRP", "TKN","TN", "TP", "DIN","TON")
cq.dat=merge(subset(q.sle.dat.xtab,WY>2009),subset(wq.dat.xtab,WY>2009&Station.ID=="C44S80")[,wq.vars],c("Date.EST","WY"))
cq.dat$log.S80.QQmean=log(cq.dat$S80.QQmean)
cq.dat$log.S80=log(cq.dat$S80+1)
cq.dat$log.TP=log(cq.dat$TP)
cq.dat$log.TN=log(cq.dat$TN)
head(cq.dat)

plot(TP~C44.QQmean,cq.dat)
plot(TN~C44.QQmean,cq.dat)

with(cq.dat,cor.test(TN,S80,method="spearman"))
with(cq.dat,cor.test(TP,S80,method="spearman"))
plot(TON~S80.QQmean,cq.dat,log="xy")
plot(DIN~S80.QQmean,cq.dat,log="xy")
plot(TON~C44.QQmean,cq.dat,log="xy")
plot(DIN~C44.QQmean,cq.dat,log="xy")


WYs=2010:2019
cols=wesanderson::wes_palette("Zissou1",length(WYs),"continuous")
#tiff(filename=paste0(plot.path,"CRE_S80CQ.tiff"),width=6.5,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.1,0.1),oma=c(3,2.5,1,0.5));
layout(matrix(c(1,2,3,3),2,2),widths=c(1,0.25))
xlim.val=c(0.01,12);xmaj=log.scale.fun(xlim.val,"major");xmin=log.scale.fun(xlim.val,"minor")

ylim.val=c(0.05,0.600);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
plot(TP~S80.QQmean,cq.dat,log="xy",type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYs)){
  with(subset(cq.dat,WY==WYs[i]),points(S80.QQmean,TP,pch=21,lwd=0.01,col=adjustcolor("black",0.5),bg=adjustcolor(cols[i],0.5),cex=1.25))
}
axis_fun(1,xmaj,xmin,NA)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=2,line=2.5,"TP (\u03BCg L\u207B\u00B9)")
mtext(side=3,"S-80")

ylim.val=c(0.5,3);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
plot(TN~S80.QQmean,cq.dat,log="xy",type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYs)){
  with(subset(cq.dat,WY==WYs[i]),points(S80.QQmean,TN,pch=21,lwd=0.01,col=adjustcolor("black",0.5),bg=adjustcolor(cols[i],0.5),cex=1.25))
}
axis_fun(1,xmaj,xmin,format(xmaj,scientific=F))
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"TN (mg L\u207B\u00B9)")
mtext(side=1,line=2,expression(paste("Q"["S80"]," : Q"[bar("X")])))

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,WYs,
       pch=21,lty=NA,lwd=0.1,
       col=adjustcolor("black",0.5),
       pt.bg=adjustcolor(cols,0.5),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,title.adj = 0,title="Water Year")
dev.off()

# monthly conc model ------------------------------------------------------
sle.hydro.mon=merge(q.sle.dat.xtab.mon,lakeO.stg.mon,c("monCY","WY"))
sle.hydro.wq.mon=merge(sle.hydro.mon,wq.dat.xtab.mon,c("WY","monCY"))

sle.hydro.wq.mon$hydro.season=as.numeric(FL.Hydroseason(sle.hydro.wq.mon$monCY)=="B_Dry")
sle.hydro.wq.mon$month=as.numeric(format(sle.hydro.wq.mon$monCY,"%m"))
sle.hydro.wq.mon=merge(sle.hydro.wq.mon,data.frame(month=c(5:12,1:4),month.plot=1:12),"month",all.x=T)
sle.hydro.wq.mon=sle.hydro.wq.mon[order(sle.hydro.wq.mon$monCY),]
sle.hydro.wq.mon$mean.TP.ugL=sle.hydro.wq.mon$mean.TP*1000

#sle.hydro.wq.mon=subset(sle.hydro.wq.mon,WY>2009)
sle.hydro.wq.mon=subset(sle.hydro.wq.mon,WY>2011) #modified due to data gap between 2010 and 2012
plot(mean.TP~monCY,subset(sle.hydro.wq.mon,WY>2011))

# data explore ------------------------------------------------------------

#tiff(filename=paste0(plot.path,"sle_scatterplot_month.tiff"),width=7,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.1,0.1),oma=c(3,3.5,0.75,0.5));
layout(matrix(1:49,7,7))

params=c("mean.TP.ugL","mean.TN","S308","S80","C44Basin", "basin.q.ratio","LakeO","del.stg")
summary(sle.hydro.wq.mon[,params])
axis.lab=c("S80 TP\n(\u03BCg L\u207B\u00B9)",
           "S80 TN\n(mg L\u207B\u00B9)",
           "Q S308\n(km\u00B3 mon\u207B\u00B9)",
           "Q S80\n(km\u00B3 mon\u207B\u00B9)",
           "Q C44\n(km\u00B3 mon\u207B\u00B9)",
           "Basin Ratio\n(Ratio)",
           "Lake O Stg\n(Ft)",
           "\u0394 Lake Stage\n(Ft)")

for(j in 1:7){
  if(j!=1){for(k in 1:(j-1)){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}}
  
  params2=params[-1:-j]
  axis.lab2=axis.lab[-1:-j]
  lim.min=c(0,0,0,0,0,0,9.5,-0.1)
  lim.max=c(500,3,0.25,0.40,0.30,1,17,0.1);by.val=c(200,1,0.1,0.2,0.1,0.5,2,0.05)
  for(i in 1:length(params2)){
    xlim.val=c(lim.min[j],lim.max[j]);by.x=by.val[j];xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
    ylim.val=c(lim.min[-1:-j][i],lim.max[-1:-j][i]);by.y=by.val[-1:-j][i];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
    xmaj.lab=if(max(xmaj)>1e4){xmaj/1000}else{xmaj}
    ymaj.lab=if(max(ymaj)>1e4){ymaj/1000}else{ymaj}
    plot(sle.hydro.wq.mon[,params[j]],sle.hydro.wq.mon[,params2[i]],xlim=xlim.val,ylim=ylim.val,axes=F,type="n",ylab=NA,xlab=NA)
    abline(h=ymaj,v=xmaj,lty=3,col="grey")
    points(sle.hydro.wq.mon[,params[j]],sle.hydro.wq.mon[,params2[i]],pch=21,bg=adjustcolor("dodgerblue1",0.25),col=adjustcolor("grey",0.5),lwd=0.01,cex=1)
    y.val=na.omit(sle.hydro.wq.mon[,c(params[j],params2[i])])[,params2[i]]
    x.val=na.omit(sle.hydro.wq.mon[,c(params[j],params2[i])])[,params[j]]
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

cor.mon=rcorr(as.matrix(sle.hydro.wq.mon[,params]),type="spearman")
cor.mon$r=round(cor.mon$r,2)
cor.mon$r[upper.tri(cor.mon$r,diag=F)]=NA
cor.mon$P=with(cor.mon,ifelse(P<0.05,"<0.05",round(P,2)))
cor.mon$P[upper.tri(cor.mon$P,diag=F)]=NA
cor.mon2=matrix(with(cor.mon,paste0(r," (",P,")")),ncol=length(params))
diag(cor.mon2)=NA
cor.mon2[cor.mon2=="NA (NA)"]=NA
param.names=c("TP","TN","Q S308","Q S80","Q C44","Basin Ratio","Lake O", "\u0394 Lake Stage")
rownames(cor.mon2)=param.names
colnames(cor.mon2)=param.names

knitr::kable(cor.mon2[2:8,1:7],align=c("c"),escape=F)%>%
  kable_styling( full_width = F)

# No PCA

# Hydrodynamic WQ model ---------------------------------------------------

# TP Model ----------------------------------------------------------------
vars=c("mean.TP","C44Basin","S308","S80","basin.q.ratio","LakeO","del.stg")
S80.TP.mod.all=lm(log(mean.TP)~.,na.omit(sle.hydro.wq.mon[,vars]))
#AIC model
S80.TP.mod.sw=stepAIC(S80.TP.mod.all,direction="both",trace=F)
S80.TP.mod.sw$anova
summary(S80.TP.mod.sw)
formula(S80.TP.mod.sw)
vif(S80.TP.mod.sw)

S80.TP.mod=lm(log(mean.TP)~C44Basin+S308+LakeO+del.stg,sle.hydro.wq.mon)
layout(matrix(1:4,2,2));plot(S80.TP.mod)
shapiro.test(residuals(S80.TP.mod));hist(residuals(S80.TP.mod))
acf(S80.TP.mod$residuals)
vif(S80.TP.mod)
summary(S80.TP.mod)

## Train vs test
# http://r-statistics.co/Linear-Regression.html
# using random 80% of data
# range(sle.hydro.wq.mon2$WY)
# sle.hydro.wq.mon2=subset(sle.hydro.wq.mon2,WY%in%seq(2010,2018,1))
set.seed(1)
tr.index=sample(1:nrow(sle.hydro.wq.mon),nrow(sle.hydro.wq.mon)*0.7)

sle.hydro.wq.mon[tr.index,"index"]="train"
sle.hydro.wq.mon[-tr.index,"index"]="test"
boxplot(mean.TP~index,sle.hydro.wq.mon)
kruskal.test(mean.TP~index,sle.hydro.wq.mon)

plot(mean.TP~monCY,sle.hydro.wq.mon)
with(sle.hydro.wq.mon[tr.index,],points(monCY,mean.TP,pch=21,bg="red"))

mod.TP=lm(log(mean.TP)~C44Basin+LakeO+del.stg,sle.hydro.wq.mon[tr.index,])
layout(matrix(1:4,2,2));plot(mod.TP)
shapiro.test(residuals(mod.TP))
acf(mod.TP$residuals)
vif(mod.TP)
summary(mod.TP)

gvlma::gvlma(mod.TP)
lmtest::bgtest(mod.TP)

mod.TP.pred=predict(mod.TP,sle.hydro.wq.mon[-tr.index,],interval="confidence")

actuals_preds <-data.frame(cbind(actuals=sle.hydro.wq.mon[-tr.index,"mean.TP"],predicted=exp(mod.TP.pred[,1])))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#tiff(filename=paste0(plot.path,"C44TP_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(0.075,0.4);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
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

#min_max_accuracy
mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
#rcompanion::accuracy(list(mod.TP))

dev.off()
plot(mean.TP~monCY,sle.hydro.wq.mon)
with(sle.hydro.wq.mon,pt_line(monCY,mean.TP,2,"grey50",1,21,"dodgerblue1",pt.col="grey50",pt.lwd=0.1))
mod.pred2=predict(mod.TP,sle.hydro.wq.mon,interval="confidence")
shaded.range(sle.hydro.wq.mon$monCY,exp(mod.pred2[,2]),exp(mod.pred2[,3]),"indianred1",lty=1)
lines(sle.hydro.wq.mon$monCY,exp(mod.pred2[,1]),col="indianred1",lwd=1.5)

#k-folding
#http://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/
library(caret)
kmodel=train(
  formula(mod.TP),
  sle.hydro.wq.mon,
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

kfold.mod.pred=predict(kmodel,na.omit(sle.hydro.wq.mon),interval="confidence")
lines(na.omit(sle.hydro.wq.mon)$monCY,exp(kfold.mod.pred),col="black",lwd=2)

kmodel2=data.frame()
kerrors=data.frame()
k=10
set.seed(1)
for(i in 1:k){
  tr.index=sample(1:nrow(sle.hydro.wq.mon),nrow(sle.hydro.wq.mon)*0.7)
  
  train.dat=sle.hydro.wq.mon[tr.index,]
  test.dat=sle.hydro.wq.mon[-tr.index,]
  
  k.cv.mod=lm(formula(mod.TP),train.dat)
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
#tiff(filename=paste0(plot.path,"C44TP_kmodel.tiff"),width=4,height=2.25,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.25,0.25));
xlim.val=c(0.075,0.4);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
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

# TN Model ----------------------------------------------------------------
plot(mean.TN~monCY,sle.hydro.wq.mon)
vars=c("mean.TN","C44Basin","S308","S80","basin.q.ratio","LakeO","del.stg")
S80.TN.mod.all=lm(log(mean.TN)~.,na.omit(sle.hydro.wq.mon[,vars]))
layout(matrix(1:4,2,2));plot(S80.TN.mod.all)
shapiro.test(residuals(S80.TN.mod.all))

dev.off()
#AIC model
S80.TN.mod.sw=stepAIC(S80.TN.mod.all,direction="both",trace=F)
S80.TN.mod.sw$anova
summary(S80.TN.mod.sw)
vif(S80.TN.mod.sw)

## Train vs test
set.seed(1)
tr.index=sample(1:nrow(sle.hydro.wq.mon),nrow(sle.hydro.wq.mon)*0.7)

sle.hydro.wq.mon[tr.index,"index"]="train"
sle.hydro.wq.mon[-tr.index,"index"]="test"

boxplot(mean.TN~index,sle.hydro.wq.mon)
kruskal.test(mean.TN~index,sle.hydro.wq.mon)

plot(mean.TN~monCY,sle.hydro.wq.mon)
with(sle.hydro.wq.mon[tr.index,],points(monCY,mean.TN,pch=21,bg="red"))

dev.off()

mod.TN=lm(log(mean.TN) ~ S308+C44Basin+LakeO ,sle.hydro.wq.mon[tr.index,])
summary(mod.TN)
layout(matrix(1:4,2,2));plot(mod.TN)
shapiro.test(residuals(mod.TN))
acf(mod.TN$residuals)
vif(mod.TN)

gvlma::gvlma(mod.TN)
lmtest::bgtest(mod.TN)

mod.TN.pred=predict(mod.TN,sle.hydro.wq.mon[-tr.index,],interval="confidence")

actuals_preds <-data.frame(cbind(actuals=sle.hydro.wq.mon[-tr.index,"mean.TN"],predicted=exp(mod.TN.pred[,1])))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#tiff(filename=paste0(plot.path,"C44TN_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(1,2);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
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
mtext(side=1,line=1.75,"Predicted TP (\u03BCg L\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TP (\u03BCg L\u207B\u00B9)")
legend("topleft",legend=c("Training Dataset","1:1 line","Correlation"),
       pch=c(21,NA,NA),lty=c(NA,1,1),lwd=c(0.1,1,2),
       col=c(adjustcolor("black",0.5),"black","indianred1"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()
