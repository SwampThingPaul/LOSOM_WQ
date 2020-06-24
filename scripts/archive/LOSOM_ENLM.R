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
#stg.dbkeys=data.frame(SITE=c("S79_H","S79_H","S235_T","S235_T","S77_H","S77_H"),DBKEY=c("00864","AN786","15566","38259","J8188","00852"))
#stg.dbkeys=subset(stg.dbkeys,DBKEY!="00864");#not sure of the datum (NGVD vs NAVD)
stg.dbkeys=data.frame(SITE=c("S79_H","S235_T","S235_T","S77_T"),DBKEY=c("AN786","15566","38259","JI497"))
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
stg.cre.dat$Date.EST=date.fun(stg.cre.dat$Date)
range(subset(stg.cre.dat,SITE=="S77_T")$Data.Value)

stg.dbkeys.bk=data.frame(SITE="S77_T",DBKEY="AI539")
stg.dat.bk=DBHYDRO_breakpoint(dates[1],dates[2],stg.dbkeys.bk$DBKEY)
stg.dat.bk=merge(stg.dat.bk,stg.dbkeys.bk,"DBKEY")
stg.dat.bk$Date.EST=date.fun(stg.dat.bk$DATE)
range(stg.dat.bk$Data.Value)
range(subset(stg.dat.bk,Data.Value>0)$Data.Value)

stg.cre.dat2=rbind(stg.cre.dat[,c("SITE","Date.EST","Data.Value")],
                   ddply(subset(stg.dat.bk,Data.Value>8),c("SITE","Date.EST"),summarise,Data.Value=mean(Data.Value,na.rm=T)))

stg.dat.da.xtab=cast(stg.cre.dat2,Date.EST~SITE,value="Data.Value",mean)
plot(S77_T~S235_T,stg.dat.da.xtab);abline(0,1,col="red")

stg.dat.da.xtab$WY=WY(stg.dat.da.xtab$Date)
stg.dat.da.xtab$month=format(stg.dat.da.xtab$Date,"%m")
stg.dat.da.xtab$CY=format(stg.dat.da.xtab$Date,"%Y")
stg.dat.da.xtab$monCY=with(stg.dat.da.xtab,date.fun(paste(CY,month,"01",sep="-")))
stg.dat.da.xtab$grad=with(stg.dat.da.xtab,(S77_T-S79_H))
stg.dat.da.xtab$grad2=with(stg.dat.da.xtab,(S235_T-S79_H))
#stg.dat.da.xtab$grad2=with(stg.dat.da.xtab,(S77_H-S79_H)); # Limited daily stage POR
#stg.dat.da.xtab$grad=with(stg.dat.da.xtab,(S235_T-S79_H)/213356.3)

plot(grad~Date.EST,stg.dat.da.xtab)
with(stg.dat.da.xtab,lines(Date.EST,grad2,col="red"))

# Sanity Check
# head(stg.dat.da.xtab)
# plot(grad2~Date,stg.dat.da.xtab,type="l")
# with(stg.dat.da.xtab,lines(Date,grad,col="red"))

# Monthly aggregation
stg.dat.da.xtab.mon=ddply(stg.dat.da.xtab,c("monCY","WY"),summarise,
                          S79_H=mean(S79_H,na.rm=T),
                          S77_T=mean(S77_T,na.rm=T),
                          S235_T=mean(S235_T,na.rm=T),
                          grad=mean(grad,na.rm=T))


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
wq.dat.xtab.mon=ddply(subset(wq.dat.xtab,Station.ID=="S79"&abs(Q.cfs)>0),c("WY","monCY"),summarise,mean.TP=mean(TP,na.rm=T),mean.TN=mean(TN,na.rm=T),mean.DIN=mean(DIN,na.rm=T),mean.TON=mean(TON,na.rm=T))
wq.dat.xtab.mon
# 
# # Doering Model
# WYs=seq(2010,2014,1)
# #WYs=seq(1999,2019,1)
# wq.vars=c("Station.ID","WY", "Date.EST", "Chla", "NH4", "NOx", "SRP", "TKN","TN", "TP", "DIN")
# doering.dat1=merge(subset(q.cre.dat.xtab,WY%in%WYs),subset(wq.dat.xtab,WY%in%WYs&Station.ID=="S79")[,wq.vars],c("Date.EST","WY"))
# head(doering.dat1)
# 
# doering.dat1$TP_kgd=with(doering.dat1,Load.Calc.kg(S79,TP))
# #doering.dat1$TP_kgd.detrend=doering.dat1$TP_kgd-mean(doering.dat1$TP_kgd,na.rm=T)
# doering.dat1$TN_kgd=with(doering.dat1,Load.Calc.kg(S79,TN))
# 
# ## TP Model
# doering.mod.TP=lm(TP_kgd~C43.basin.in+LakeQ,doering.dat1)
# 
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


# monthly conc model ------------------------------------------------------
cre.hydro.mon=merge(q.cre.dat.xtab.mon,stg.dat.da.xtab.mon,c("monCY","WY"))
cre.hydro.wq.mon=merge(cre.hydro.mon,wq.dat.xtab.mon,c("WY","monCY"))

cre.hydro.wq.mon$hydro.season=as.numeric(FL.Hydroseason(cre.hydro.wq.mon$monCY)=="B_Dry")
cre.hydro.wq.mon$month=as.numeric(format(cre.hydro.wq.mon$monCY,"%m"))
cre.hydro.wq.mon=merge(cre.hydro.wq.mon,data.frame(month=c(5:12,1:4),month.plot=1:12),"month",all.x=T)
cre.hydro.wq.mon=cre.hydro.wq.mon[order(cre.hydro.wq.mon$monCY),]
cre.hydro.wq.mon$mean.TP.ugL=cre.hydro.wq.mon$mean.TP*1000

cre.hydro.wq.mon=subset(cre.hydro.wq.mon,WY>2009)

vars=c("mean.TP","grad","C43Basin","S77","S78","basin.q.ratio","S79_H","S235_T","WY")
#tiff(filename=paste0(plot.path,"CRE_scatterplot_month.tiff"),width=7,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.1,0.1),oma=c(3,3.5,0.75,0.5));
layout(matrix(1:49,7,7))

params=c("mean.TP.ugL","mean.TN","S77","S79","S78Qfreq","C43Basin", "basin.q.ratio","grad")
summary(cre.hydro.wq.mon[,params])
axis.lab=c("S79 TP\n(\u03BCg L\u207B\u00B9)",
           "S79 TN\n(mg L\u207B\u00B9)",
           "Q S77\n(km\u00B3 mon\u207B\u00B9)",
           "Q S79\n(km\u00B3 mon\u207B\u00B9)",
           "Q S78>0 Freq\n(days)",
           "Q C43\n(km\u00B3 mon\u207B\u00B9)",
           "Basin Ratio\n(Ratio)",
           "Mean Grad\n(Ft)")

for(j in 1:7){
  if(j!=1){for(k in 1:(j-1)){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}}
  
  params2=params[-1:-j]
  axis.lab2=axis.lab[-1:-j]
  lim.min=c(0,0,0,0,0,0,0,5.6)
  lim.max=c(300,3,0.50,1,30,0.60,1,8.6);by.val=c(150,1.5,0.25,0.5,15,0.3,0.5,1.5)
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
param.names=c("TP","TN","Q S77","Q S79","Q S78 Freq","Q C43","Basin Ratio","Gradient")
rownames(cor.mon2)=param.names
colnames(cor.mon2)=param.names

knitr::kable(cor.mon2[2:8,1:7],align=c("c"),escape=F)%>%
  kable_styling( full_width = F)

# Quick PCA
library(REdaS)
library(vegan)
params=c("mean.TP.ugL","mean.TN","mean.DIN","mean.TON","S77","S78","S78Qfreq","C43Basin", "basin.q.ratio","grad")
my.pca=rda(na.omit(cre.hydro.wq.mon[,params]),scale=T)

plot(my.pca)

eig <- my.pca$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)

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
           "TN","DIN","TON",
           "Q S77",
           "Q S78",
           "Q S78>0 Freq",
           "Q C43",
           "Basin Ratio",
           "Stg Grad")

#tiff(filename=paste0(plot.path,"CRE_month_PCA_biplot.tiff"),width=4.5,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,3,0.1,0.5),oma=c(1.5,1.5,0.75,0.5));
#layout(matrix(1:2,2,1))

xlim.val=c(-1,2.5);by.x=1;xmaj=c(seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-2.5,2.5);by.y=1.25;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs$sites[,c(1,2)],pch=21,bg=adjustcolor("grey",0.20),col=adjustcolor("black",0.20),cex=1,lwd=0.1); #plots the points
arrows(0,0,scrs$species[,1],scrs$species[,2],length = 0.07, angle = 25, code = 2,col="indianred1",lwd=2);# makes the arrows
with(scrs,text(species[,1],species[,2],labels=axis.lab,cex=1,font=2,pos=1));#adds labels to the arrows; 
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

# S79 TP concentration ----------------------------------------------------

plot(mean.TP~monCY,cre.hydro.wq.mon,col="red")

# Time series stationarity
#https://rpubs.com/richkt/269797
plot(mean.TP~monCY,cre.hydro.wq.mon,type="l",col="red")
acf(cre.hydro.wq.mon$mean.TP)

Box.test(cre.hydro.wq.mon$mean.TP,type="Ljung-Box")#non-stationary p<0.05
tseries::adf.test(cre.hydro.wq.mon$mean.TP)#no unit-root p<0.05
tseries::kpss.test(cre.hydro.wq.mon$mean.TP,null="Trend")# non-stationary and no unit-root p>0.05

cre.hydro.wq.mon$decWY=decimal.WY(cre.hydro.wq.mon$monCY)
TP.trend=lm(mean.TP~decWY,cre.hydro.wq.mon)
plot(TP.trend$residuals,type="o")
acf(TP.trend$residuals)

plot(diff(cre.hydro.wq.mon$mean.TP))
acf(diff(cre.hydro.wq.mon$mean.TP))


#seasonal analysis
ylim.val=c(0.05,0.35);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.3);ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1 - 0.1, 12 + 0.4 + 0.5);xmaj=1:12;xmaj.lab=c(5:12,1:4)
WYs=seq(2009,2019,1)
cols=wesanderson::wes_palette("Zissou1",length(WYs),"continuous")
#tiff(filename=paste0(plot.path,"S79_monthlyTP_years.tiff"),width=6.5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,0.75,0.25,0.1),oma=c(3,2.75,1,0.5));
plot(mean.TP~month.plot,cre.hydro.wq.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYs)){
with(subset(cre.hydro.wq.mon,WY==WYs[i]),pt_line(month.plot,mean.TP,2,cols[i],1,21,cols[i],pt.lwd=0.1))
text(12,subset(cre.hydro.wq.mon,WY==WYs[i]&month.plot==12)$mean.TP,WYs[i],col=cols[i],pos=4)
}
axis_fun(1,line=-0.5,xmaj,xmaj,xmaj.lab)
axis_fun(2,ymaj,ymin,format(ymaj*1000));box(lwd=1)
mtext(side=3,"S-79")
mtext(side=1,line=1.5,"Month")
mtext(side=2,line=2.5,"Total Phosphorus (\u03BCg L\u207B\u00B9)")
dev.off()

cre.hydro.wq.mon$CY=as.numeric(format(cre.hydro.wq.mon$monCY,"%Y"))
#tiff(filename=paste0(plot.path,"S79_TPseasonal.tiff"),width=7,height=1.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,0.75,0.25,0.1),oma=c(2,2.75,1,0.5));
layout(matrix(1:12,1,12,byrow=T))

ylim.val=c(0.05,0.30);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.3);ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(2008,2020);by.x=4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)

for(i in c(5:12,1:4)){
  tmp.dat=subset(cre.hydro.wq.mon,month==i)
  
  plot(mean.TP~CY,cre.hydro.wq.mon,xlim=xlim.val,ylim=ylim.val,axes=F,ylab=NA,xlab=NA,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(tmp.dat,points(CY,mean.TP,pch=21,bg="grey",lwd=0.01,cex=0.75))
  mod=mblm(mean.TP~CY,na.omit(tmp.dat[,c("mean.TP","CY")]),repeated = F)
  x.val=seq(min(tmp.dat$CY),max(tmp.dat$CY),length.out=20)
  mod.pred=predict(mod,data.frame(CY=x.val),interval="confidence")
  lines(x.val,mod.pred[,1],lwd=1,col="red")
  lines(x.val,mod.pred[,2],lty=2,col="red")
  lines(x.val,mod.pred[,3],lty=2,col="red")
  if(i==5){axis_fun(2,ymaj,ymin,ymaj*1000,cex=0.8,maj.tcl=-0.3,min.tcl=-0.15,lwd=0.8,line=-0.5)}else{axis_fun(2,ymaj,ymin,NA,cex=0.8,maj.tcl=-0.3,min.tcl=-0.15,lwd=0.8)}
  axis_fun(1,line=-1,xmaj,xmin,xmaj,cex=0.6,maj.tcl=-0.3,min.tcl=-0.15,lwd=0.8)
  box(lwd=0.8)
  mtext(side=3,month.abb[i],cex=0.75)
}  
mtext(side=1,line=0.7,outer=T,"Calendar Year",cex=0.6)
mtext(side=2,line=1.25,outer=T,"Total Phosphorus (\u03BCg L\u207B\u00B9)",cex=0.6)
dev.off()

# range(cre.hydro.wq.mon$monCY)
# TP.ts=ts(cre.hydro.wq.mon$mean.TP,start=c(2009,05),frequency = 12)
# decomp.TP.ts=decompose(TP.ts,type="mult")
# plot(decomp.TP.ts)
# library(forecast)
# decomp.TP.ts2=stl(TP.ts,s.window="periodic")
# decomp.TP.ts2
# plot(decomp.TP.ts2$time.series)

#seasonal kendall 
# Month is season
with(cre.hydro.wq.mon,rkt(lubridate::decimal_date(monCY),mean.TP,month))
# No seasonal (monthly) trend

EnvStats::kendallSeasonalTrendTest(y=cre.hydro.wq.mon$mean.TP,season=cre.hydro.wq.mon$month,year=cre.hydro.wq.mon$CY)

#assumptions
bartlett.test(log(mean.TP)~as.factor(month),cre.hydro.wq.mon)
car::leveneTest(log(mean.TP)~as.factor(month),cre.hydro.wq.mon)
### Welch’s anova for unequal variances - do months significantly vary?
oneway.test(log(mean.TP)~as.factor(month),cre.hydro.wq.mon, var.equal = F)
#months signifiantly differ

#tiff(filename=paste0(plot.path,"S79_monthlyTP.tiff"),width=6.5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,0.75,0.25,0.1),oma=c(3,2.75,1,0.5));

ylim.val=c(0.05,0.35);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.3);ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1,12);xmaj=1:12;xmaj.lab=c(5:12,1:4)
plot(mean.TP~month.plot,cre.hydro.wq.mon,log="y",ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(month.plot,mean.TP,pch=21,bg=adjustcolor("dodgerblue1",0.25),col="dodgerblue1",lwd=0.1,cex=1.25))
k=predict(loess(mean.TP~month.plot,cre.hydro.wq.mon),data.frame(month.plot=seq(1,12,length.out = 24)),se=T)
lines(seq(1,12,length.out = 24),k$fit,lwd=2,col="red")
lines(seq(1,12,length.out = 24),k$fit - qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
lines(seq(1,12,length.out = 24),k$fit + qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
axis_fun(1,line=-0.5,xmaj,xmaj,xmaj.lab)
axis_fun(2,ymaj,ymin,format(ymaj*1000));box(lwd=1)
mtext(side=3,"S-79")
mtext(side=1,line=2,"Month")
mtext(side=2,line=2.5,"Total Phosphorus (\u03BCg L\u207B\u00B9)")
legend("topright",c("Monthly Mean","LOESS \u00B1 95% CI"),
       pch=c(21,NA),
       lty=c(NA,1),lwd=c(0.01,1),
       col=c("dodgerblue1","red"),
       pt.bg=c(adjustcolor("dodgerblue1",0.25),NA),
       pt.cex=1,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()


ylim.val=c(1,2);ymaj=c(0.05,log.scale.fun(ylim.val,"major"),0.3);ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1,12);xmaj=1:12;xmaj.lab=c(5:12,1:4)
plot(mean.TN~month.plot,cre.hydro.wq.mon,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(cre.hydro.wq.mon,points(month.plot,mean.TN,pch=21,bg=adjustcolor("dodgerblue1",0.25),col="dodgerblue1",lwd=0.1,cex=1.25))
k=predict(loess(mean.TN~month.plot,cre.hydro.wq.mon),data.frame(month.plot=seq(1,12,length.out = 24)),se=T)
lines(seq(1,12,length.out = 24),k$fit,lwd=2,col="red")
lines(seq(1,12,length.out = 24),k$fit - qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
lines(seq(1,12,length.out = 24),k$fit + qt(0.975,k$df)*k$se, lty=2,lwd=1,col="red")
axis_fun(1,line=-0.5,xmaj,xmaj,xmaj.lab)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)


# Hydrodynamic ------------------------------------------------------------
null.mod=lm(mean.TP~1,cre.hydro.wq.mon)
summary(null.mod)
layout(matrix(1:4,2,2));plot(null.mod)
shapiro.test(residuals(null.mod))
acf(residuals(null.mod))
#car::vif(null.mod)

#vars=c("mean.TP","grad","C43Basin","S78","S79","S77","basin.q.ratio","hydro.season")
#vars=c("mean.TP","grad","S79_H","S235_T","C43Basin","S77","S78","basin.q.ratio","S78Qfreq")
vars=c("mean.TP","grad","S79_H","S77_T","C43Basin","S77","S78","basin.q.ratio","S78Qfreq")
S79.TP.mod.all=lm(log(mean.TP)~.,na.omit(cre.hydro.wq.mon[,vars]))
#AIC model
S79.TP.mod.sw=stepAIC(S79.TP.mod.all,direction="both",trace=F)
S79.TP.mod.sw$anova
summary(S79.TP.mod.sw)
formula(S79.TP.mod.sw)

S79.TP.mod=lm(formula(S79.TP.mod.sw),cre.hydro.wq.mon)
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

mod=lm(formula(S79.TP.mod.sw),cre.hydro.wq.mon[tr.index,])
#mod=lm(log(mean.TP)~grad+C43Basin+S78Qfreq+basin.q.ratio+S79_H,cre.hydro.wq.mon[tr.index,])
layout(matrix(1:4,2,2));plot(mod)
shapiro.test(residuals(mod))
acf(mod$residuals)
vif(mod)
summary(mod)
AIC(mod)

gvlma::gvlma(mod)
lmtest::bgtest(mod)

broom::glance(mod)
broom::tidy(mod)

mod.pred=predict(mod,cre.hydro.wq.mon[-tr.index,],interval="confidence")

actuals_preds <-data.frame(cbind(actuals=cre.hydro.wq.mon[-tr.index,"mean.TP"],predicted=exp(mod.pred[,1])))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#https://rcompanion.org/handbook/G_14.html
#min_max_accuracy
mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  

dev.off()
plot(mean.TP~monCY,cre.hydro.wq.mon)
with(cre.hydro.wq.mon,pt_line(monCY,mean.TP,2,"grey50",1,21,"dodgerblue1",pt.col="grey50",pt.lwd=0.1))
mod.pred2=predict(mod,cre.hydro.wq.mon,interval="confidence")
shaded.range(cre.hydro.wq.mon$monCY,exp(mod.pred2[,2]),exp(mod.pred2[,3]),"indianred1",lty=1)
lines(cre.hydro.wq.mon$monCY,exp(mod.pred2[,1]),col="indianred1",lwd=1.5)
#lines($monCY,exp(mod.pred2[,2]),col="indianred1",lty=2)
#lines($monCY,exp(mod.pred2[,3]),col="indianred1",lty=2)

#k-folding
library(caret)
kmodel=train(
  formula(S79.TP.mod.sw),
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

# library(DAAG)
# cv.rslt=CVlm(data=cre.hydro.wq.mon,
#              form.lm=log(mean.TP)~grad+C43Basin+S78Qfreq+basin.q.ratio+S79_H,
#              m=5,seed=1,
#              plotit="Observed")
# CVlm(data=cre.hydro.wq.mon,
#              form.lm=log(mean.TP)~grad+C43Basin+S78Qfreq+basin.q.ratio+S79_H,
#              m=5,seed=1,
#              plotit="Residual")
# 
# attr(cv.rslt,'ms');# mean square error

## TN
plot(mean.TN~monCY,cre.hydro.wq.mon)

#subset(cre.hydro.wq.mon,mean.TN>2)
#subset(wq.dat.xtab,month=="06"&CY==2011&Station.ID=="S79")

#cre.hydro.wq.mon$mean.TN=with(cre.hydro.wq.mon,ifelse(mean.TN>2,NA,mean.TN))

#vars=c("mean.TN","grad","S79_H","S235_T","C43Basin","S77","S78","basin.q.ratio","S78Qfreq")
vars=c("mean.TN","grad","S79_H","S77_T","C43Basin","S77","S78","basin.q.ratio","S78Qfreq")
S79.TN.mod.all=lm(log(mean.TN)~.,na.omit(cre.hydro.wq.mon[,vars]))
layout(matrix(1:4,2,2));plot(S79.TN.mod.all)
shapiro.test(residuals(S79.TN.mod.all))

dev.off()
#AIC model
S79.TN.mod.sw=stepAIC(S79.TN.mod.all,direction="both",trace=F)
S79.TN.mod.sw$anova
summary(S79.TN.mod.sw)
vif(S79.TN.mod.sw)
# #StepVIF
# library(pedometrics)
# S79.TN.mod.sVIF=stepVIF(S79.TN.mod.all)
# summary(S79.TN.mod.sVIF)
# vif(S79.TN.mod.sVIF)

## Train vs test
set.seed(1)
tr.index=sample(1:nrow(cre.hydro.wq.mon),nrow(cre.hydro.wq.mon)*0.7)

cre.hydro.wq.mon[tr.index,"index"]="train"
cre.hydro.wq.mon[-tr.index,"index"]="test"
boxplot(mean.TN~index,cre.hydro.wq.mon)
kruskal.test(mean.TN~index,cre.hydro.wq.mon)

plot(mean.TN~monCY,cre.hydro.wq.mon)
with(cre.hydro.wq.mon[tr.index,],points(monCY,mean.TN,pch=21,bg="red"))

#mod.TN=lm(log(mean.TN)~S79_H+S77+S78+S78Qfreq,cre.hydro.wq.mon[tr.index,])
mod.TN=lm(formula(S79.TN.mod.sw),cre.hydro.wq.mon[tr.index,])
layout(matrix(1:4,2,2));plot(mod.TN)
shapiro.test(residuals(mod.TN))
acf(mod.TN$residuals)
vif(mod.TN)
summary(mod.TN)
AIC(mod.TN)

plot(S77~S78,cre.hydro.wq.mon[tr.index,])
plot(mean.TN~C43Basin,cre.hydro.wq.mon[tr.index,])

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


# S79 hind cast models ----------------------------------------------------
### hind cast eval
hind.hydro=subset(cre.hydro.mon,WY<2010)
hind.wq=subset(wq.dat.xtab.mon,WY<2010)
hind.wq.hydro=merge(hind.wq,hind.hydro,c("WY","monCY"))

hind.TN=predict(mod.TN,hind.wq.hydro,interval="confidence")
plot(mean.TN~monCY,hind.wq.hydro)
points(hind.wq.hydro$monCY,exp(hind.TN[,1]),pch=21,bg="dodgerblue1")
mean(apply(cbind(hind.wq$mean.TN,exp(hind.TN[,1])), 1, min,na.rm=T) / apply(cbind(hind.wq$mean.TN,exp(hind.TN[,1])), 1, max,na.rm=T)) 
mean(abs((exp(hind.TN[,1]) - hind.wq$mean.TN))/hind.wq$mean.TN,na.rm=T)  

hind.TP=predict(mod,hind.wq.hydro,interval="confidence")
plot(mean.TP~monCY,hind.wq.hydro,type="b")
shaded.range(hind.wq.hydro$monCY,exp(hind.TP[,2]),exp(hind.TP[,3]),"grey")
points(hind.wq.hydro$monCY,exp(hind.TP[,1]),pch=21,bg="dodgerblue1")
mean(apply(cbind(hind.wq$mean.TP,exp(hind.TP[,1])), 1, min,na.rm=T) / apply(cbind(hind.wq$mean.TP,exp(hind.TP[,1])), 1, max,na.rm=T)) 
mean(abs((exp(hind.TP[,1]) - hind.wq$mean.TP))/hind.wq$mean.TP,na.rm=T)  



# STL ---------------------------------------------------------------------
STL.wmd.site=data.frame(SITE=c("S308TW","S80HW"),DBKEY=c("TA370","TA256"))#limited period of record.
STL.wmd.site.bk=data.frame(SITE=c("S308TW","S80HW"),DBKEY=c("AI533","AI560"))
usgs.site=data.frame(site_no=c("02276877", "02276998"),SITE=c("S308TW","S80HW"))

STL.wmd.stg=data.frame()
for(i in 1:nrow(STL.wmd.site)){
  tmp=DBHYDRO_daily(dates[1],dates[2],STL.wmd.site$DBKEY[i])
  tmp$DBKEY=as.character(STL.wmd.site$DBKEY[i])
  STL.wmd.stg=rbind(STL.wmd.stg,tmp)
  print(i)
}
STL.wmd.stg$Date.EST=date.fun(STL.wmd.stg$Date)
STL.wmd.stg=merge(STL.wmd.stg,STL.wmd.site,"DBKEY")
unique(STL.wmd.stg$Station)
#S308TW.stg=ddply(subset(STL.wmd.stg,SITE=="S308TW"),"Date.EST",summarise,WMD_S308TW=mean(Data.Value,na.rm=T))

STL.wmd.stg.bk=data.frame()
for(i in 1:nrow(STL.wmd.site.bk)){
  tmp=DBHYDRO_breakpoint(dates[1],dates[2],STL.wmd.site.bk$DBKEY[i])
  tmp$DBKEY=as.character(STL.wmd.site.bk$DBKEY[i])
  STL.wmd.stg.bk=rbind(STL.wmd.stg.bk,tmp)
  print(i)
}
STL.wmd.stg.bk=merge(STL.wmd.stg.bk,STL.wmd.site.bk,"DBKEY")
STL.wmd.stg.bk$Date.EST=date.fun(STL.wmd.stg.bk$DATE)
range(STL.wmd.stg.bk$Data.Value,na.rm=T)
range(subset(STL.wmd.stg.bk,Data.Value>0)$Data.Value,na.rm=T)
STL.wmd.stg.bk.da=ddply(subset(STL.wmd.stg.bk,Data.Value>0),c("SITE","Date.EST"),summarise,Data.Value=mean(Data.Value,na.rm=T))


STL.wmd.stg=rbind(STL.wmd.stg[,c("SITE","Date.EST","Data.Value")],STL.wmd.stg.bk.da)
STL.wmd.stg.xtab=data.frame(cast(STL.wmd.stg,Date.EST~SITE,value="Data.Value",mean))
STL.wmd.stg.xtab.comp=STL.wmd.stg.xtab

STL.wmd.stg.xtab$WY=WY(STL.wmd.stg.xtab$Date.EST)
STL.wmd.stg.xtab$month=format(STL.wmd.stg.xtab$Date.EST,"%m")
STL.wmd.stg.xtab$CY=format(STL.wmd.stg.xtab$Date.EST,"%Y")
STL.wmd.stg.xtab$monCY=with(STL.wmd.stg.xtab,date.fun(paste(CY,month,"01",sep="-")))
STL.wmd.stg.xtab$grad=with(STL.wmd.stg.xtab,(S308TW-S80HW))


colnames(STL.wmd.stg.xtab.comp)=c("Date.EST",paste("WMD",STL.wmd.site$SITE,sep="_"))

plot(S80HW~Date.EST,STL.wmd.stg.xtab)

##USGS compare
STL.usgs.stg=data.frame()
for(i in 1:nrow(usgs.site)){
  tmp=readNWISdv(usgs.site$site_no[i],c("00060","00065"),dates[1],dates[2])
  STL.usgs.stg=rbind(tmp,STL.usgs.stg)
}
unique(STL.usgs.stg$site_no)
STL.usgs.stg=renameNWISColumns(STL.usgs.stg)
STL.usgs.stg$Date.EST=date.fun(STL.usgs.stg$Date)
STL.usgs.stg=merge(STL.usgs.stg,usgs.site,"site_no")
STL.usgs.stg.xtab=data.frame(cast(STL.usgs.stg,Date.EST~SITE,value="GH",mean))
colnames(STL.usgs.stg.xtab)=c("Date.EST",paste("USGS",usgs.site$SITE,sep="_"))
plot(USGS_S80HW~Date.EST,STL.usgs.stg.xtab)

stg.compare=merge(STL.wmd.stg.xtab.comp,STL.usgs.stg.xtab,"Date.EST",all.y=T)

plot(WMD_S308TW~USGS_S308TW,stg.compare);abline(0,1,col="red")
plot(WMD_S80HW~USGS_S80HW,stg.compare);abline(0,1,col="red")
#double check datum conversion https://www.ngs.noaa.gov/cgi-bin/VERTCON/vert_con.prl
# for S80HW, NAVD88 - NGVD29 = -0.444 meter
summary(lm(USGS_S80HW~WMD_S80HW,stg.compare))

plot(USGS_S80HW+1.44~WMD_S80HW,stg.compare);abline(0,1,col="red")
summary(lm(USGS_S80HW+1.44~WMD_S80HW,stg.compare))
