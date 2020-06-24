## 
## LOSOM Water Quality Subteam
## Estuary Nutrient Load Model
## Caloosahatchee River Estuary (CRE)
## Dynamic model development 
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

## BAD 
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

## Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape)

library(dataRetrieval)
library(RCurl)

library(HURDAT)

# GIS libraries 
library(rgdal)
library(rgeos)
library(raster)

## Paths
wd="C:/Julian_LaCie/_Github/LOSOM_WQ"

paths=paste0(wd,c("/Plots/","/Export/","/Data/"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

GIS.path="C:/Julian_LaCie/_GISData"

## 

nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##
dates=as.Date(c("1979-05-01","2019-04-30"))

# Discharge  --------------------------------------------------------------
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
q.cre.dat$WY=WY(q.cre.dat$Date.EST)

# q.cre.dat$delta.day=with(q.cre.dat,c(NA,diff(Date.EST)))
# range(subset(q.cre.dat,SITE=="S79")$delta.day,na.rm = T)
# range(subset(q.cre.dat,SITE=="S77")$delta.day,na.rm = T)
# subset(q.cre.dat,SITE=="S79"&delta.day>1)
# subset(q.cre.dat,SITE=="S77"&delta.day!=1)

q.cre.dat.xtab=cast(q.cre.dat,Date.EST~SITE,value="Data.Value",mean)
q.cre.dat.xtab$month=format(q.cre.dat.xtab$Date,"%m")
q.cre.dat.xtab$CY=format(q.cre.dat.xtab$Date,"%Y")
q.cre.dat.xtab$monCY=with(q.cre.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
q.cre.dat.xtab$WY=WY(q.cre.dat.xtab$Date.EST)

# Sanity Check
# plot(S79~Date.EST,q.cre.dat.xtab,type="l")
# plot(S78~Date.EST,q.cre.dat.xtab,type="l")
# plot(S77~Date.EST,q.cre.dat.xtab,type="l")

#Positive discharge only
q.cre.dat.xtab$S77=with(q.cre.dat.xtab,ifelse(S77<0,0,S77))
q.cre.dat.xtab$S79=with(q.cre.dat.xtab,ifelse(S79<0,0,S79))
q.cre.dat.xtab$C43=with(q.cre.dat.xtab,ifelse(S79<S77,0,S79-S77))
# plot(C43~Date.EST,q.cre.dat.xtab,type="l")

q.cre.dat.xtab.WY=ddply(q.cre.dat.xtab,c("WY"),summarise,Q.S77=sum(cfs.to.acftd(S77),na.rm=T),Q.S79=sum(cfs.to.acftd(S79),na.rm=T))
q.cre.dat.xtab.WY$Q.C43=with(q.cre.dat.xtab.WY,ifelse(Q.S79<Q.S77,0,Q.S79-Q.S77))
q.cre.dat.xtab.WY$basin.ratio=with(q.cre.dat.xtab.WY,Q.C43/Q.S79)
# CRE WQ ------------------------------------------------------------------
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7),param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp"))
wq.sites=c("S79","S77")
wq.dat=DBHYDRO_WQ(dates[1],dates[2],wq.sites,params$Test.Number)
wq.dat=merge(wq.dat,params,"Test.Number")
unique(wq.dat$Collection.Method)
wq.dat=subset(wq.dat,Collection.Method=="G")

# plot(HalfMDL~Date.EST,subset(wq.dat,Station.ID=="S79"& param=="TP"))
# plot(HalfMDL~Date.EST,subset(wq.dat,Station.ID=="S77"& param=="TP"))

wq.dat.xtab=cast(wq.dat,Station.ID+Date.EST~param,value="HalfMDL",mean)
wq.dat.xtab$DIN=with(wq.dat.xtab,NH4+NOx)
wq.dat.xtab$TN=with(wq.dat.xtab, TN_Combine(NOx,TKN,TN))
wq.dat.xtab$TON=with(wq.dat.xtab,TN-DIN)
wq.dat.xtab$WY=WY(wq.dat.xtab$Date.EST)
wq.dat.xtab$month=format(wq.dat.xtab$Date.EST,"%m")
wq.dat.xtab$CY=format(wq.dat.xtab$Date.EST,"%Y")
wq.dat.xtab$monCY=with(wq.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
wq.dat.xtab$TNTP=with(wq.dat.xtab, (TN/14.00)/(TP/30.97))

# Reversal Evaluation
wq.dat.xtab$TPReversal=with(wq.dat.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
wq.dat.xtab$TNReversal=with(wq.dat.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

sum(wq.dat.xtab$TNReversal,na.rm=T)
sum(wq.dat.xtab$TPReversal,na.rm=T)

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,wq.dat.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,wq.dat.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

dev.off()

# Rainfall ----------------------------------------------------------------
C43Basin_rainfall.sites=data.frame(DBKEY=c(c("06093","15786"),c("VN424","15465"),c("06221","16625","06243","15495"),"06083",c("05953","06079","64082")),
                                   SITE=c(rep("PALMDALE",2),rep("WHIDDEN",2),rep("S78R",4),"KERI_TOW",rep("DEVILS",3)))
C43Basin_rainfall.sites$BASIN="C43"

rf.dat=data.frame()
pb=txtProgressBar(1,nrow(C43Basin_rainfall.sites),1,style=3)
for(i in 1:nrow(C43Basin_rainfall.sites)){
  tmp=DBHYDRO_daily(dates[1],dates[2],C43Basin_rainfall.sites$DBKEY[i])
  tmp$DBKEY=as.character(C43Basin_rainfall.sites$DBKEY[i])
  rf.dat=rbind(tmp,rf.dat)
  setTxtProgressBar(pb, i)
}
rf.dat=merge(rf.dat,C43Basin_rainfall.sites,"DBKEY")
rf.dat$Date.EST=date.fun(rf.dat$Date)
rf.dat.da=ddply(rf.dat,c("SITE","BASIN","Date.EST"),summarise,mean.RF=mean(in.to.cm(Data.Value),na.rm=T))
rf.dat.da2=ddply(rf.dat.da,c("BASIN","Date.EST"),summarise,mean.RF=mean(mean.RF,na.rm=T),N.val=N.obs(mean.RF))
rf.dat.da2$WY=WY(rf.dat.da2$Date.EST)
# plot(mean.RF~Date.EST,rf.dat.da2)
rf.dat.WY=ddply(rf.dat.da2,c("BASIN","WY"),summarise,TRF.cm=sum(mean.RF,na.rm=T))


# Lake O Stage ------------------------------------------------------------
#lakeo.stg.dbkey=data.frame(SITE=c("LakeO","LakeO"),DBKEY=c("15611","06832"),priority=c("P1","P2"));#15611 Pref DBKEY
#lakeO.stg=DBHYDRO_daily(dates[1],dates[2],lakeo.stg.dbkey$DBKEY)
lakeO.stg=DBHYDRO_daily(dates[1],dates[2],"00268")
plot(Data.Value~Date,lakeO.stg)

lakeO.stg$Date.EST=date.fun(lakeO.stg$Date)
lakeO.stg$WY=WY(lakeO.stg$Date)
#lakeO.stg$delta.stg=with(lakeO.stg,c(NA,diff(Data.Value)))
lakeO.stg.WY=ddply(lakeO.stg,"WY",summarise,mean.stg=mean(Data.Value,na.rm=T),q1.stg=quantile(Data.Value,na.rm=T,probs=0.25),q3.stg=quantile(Data.Value,na.rm=T,probs=0.75),min.stg=min(Data.Value,na.rm=T))
plot(mean.stg~WY,lakeO.stg.WY)
plot(q1.stg~WY,lakeO.stg.WY)
plot(min.stg~WY,lakeO.stg.WY)


# Hurricane analysis ------------------------------------------------------
ogrListLayers(paste0(GIS.path,"/AHED_release/AHED_20171102.gdb"))

basins=readOGR(paste0(GIS.path,"/AHED_release/AHED_20171102.gdb"),"WATERSHED")
plot(basins)

library(tmap)
tmap_mode("view")
tm_shape(basins)+tm_polygons()

CRE.basins=c("WEST CALOOSAHATCHEE","EAST CALOOSAHATCHEE",
             "NICODEMUS SLOUGH NORTH","NICODEMUS SLOUGH SOUTH","HICPOHEE NORTH",
             "S-4","LAKE OKEECHOBEE")
plot(subset(basins,NAME%in%CRE.basins))

# Hurricane Track Data
hur.dat=get_hurdat(basin = c("AL"))
range(hur.dat$DateTime)
hur.dat$date=date.fun(hur.dat$DateTime,tz="UTC")
hur.dat$Year=as.numeric(format(hur.dat$date,"%Y"))
hur.dat$WY=WY(hur.dat$DateTime)
hur.year=ddply(hur.dat,c("Key","Name"),summarise,Year.start=min(Year,na.rm=T),WY.start=min(WY,na.rm=T),max.wind.ms=max(Wind*0.44704,na.rm=T),min.pres=min(Pressure,na.rm=T))
hur.dat.sp.pt=SpatialPointsDataFrame(coords=hur.dat[,c("Lon","Lat")],data=hur.dat,proj4string = CRS("+init=epsg:4269"))

#Convert point to line data
hur.dat2=hur.dat
hur.id=ddply(hur.dat,c("Key","Name"),summarise,N.val=length(which(Key!="NA")))
hur.dat2$dLat=with(hur.dat2,ave(Lat,Key,FUN=function(x)c(0,diff(x))))
hur.dat2$dLon=with(hur.dat2,ave(Lon,Key,FUN=function(x)c(0,diff(x))))
hur.dat2$dist.m=with(hur.dat2,sqrt((dLat^2)+(dLon^2)));#calculates distance between points

hur.id=ddply(hur.dat2,c("Key","Name"),summarise,N.val=N.obs(Key),max.dist=max(dist.m,na.rm=T))
hur.id=subset(hur.id,is.na(Key)==F)
hur.dat2=subset(hur.dat2,Key%in%subset(hur.id,N.val>1&max.dist>0)$Key);#need greater than one point and some distance between points to create a line
coordinates(hur.dat2)=c("Lon","Lat")
hur.dat2=hur.dat2[order(hur.dat2$Key,hur.dat2$DateTime),]
hur.dat2$WY=WY(hur.dat2$DateTime)
path=sp::split(hur.dat2,hur.dat2$Key)

##Convert spatial points to spatial lines for each hurricane
sp_lines=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[1]])),unique(path[[1]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
pb=txtProgressBar(1,max=length(path),style=3)
for(i in 2:length(path)){
  tmp=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
  #tmp=SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269"))
  sp_lines=rbind(sp_lines,tmp)
  setTxtProgressBar(pb,i)
}
chk=data.frame(gIsValid(sp_lines,byid=T));#checks
chk$Key=rownames(chk)
colnames(chk)=c("geo.valid","Key")
subset(chk,geo.valid=="FALSE")

hur.tck=spTransform(sp_lines,utm17)
hur.tck=merge(hur.tck,hur.year,by.x=c("Key","Name"),by.y=c("Key","Name"))

#CRE boundary and buffer
plot(subset(basins,NAME%in%CRE.basins))
cre.buffer1=gBuffer(subset(basins,NAME%in%CRE.basins),width=20*1000);
plot(cre.buffer1,add=T)
cre.buffer1=SpatialPolygonsDataFrame(cre.buffer1,data.frame(row.names = "buffer",width.km=20,area="CRE"))
cre.buffer1=spTransform(cre.buffer1,utm17)

cre.hurr=over(cre.buffer1,hur.tck,returnList = T,byid=T)[[1]];#select only hurricanes that cross the 200 km buffer. 
cre.hurr=cre.hurr[order(cre.hurr$WY.start,cre.hurr$Key),]

cre.hurr=subset(cre.hurr,WY.start>=1981&max.wind.ms>10)


# Load --------------------------------------------------------------------
head(q.cre.dat)
unique(q.cre.dat$SITE)
range(q.cre.dat$Data.Value,na.rm=T)
q.cre.dat$Data.Value=with(q.cre.dat,ifelse(Data.Value<0,NA,Data.Value))
head(wq.dat.xtab)
unique(wq.dat.xtab$Station.ID)
head(subset(wq.dat.xtab,Station.ID=="S79"))
head(subset(wq.dat.xtab,Station.ID=="S77"))

q.vars=c("SITE","Date.EST","WY","Data.Value")
wq.vars=c("Station.ID","Date.EST","TP","TN")
q.wq=merge(subset(q.cre.dat[,q.vars],SITE!="S78"&WY>1981),
           wq.dat.xtab[,wq.vars],by.x=q.vars[1:2],by.y=wq.vars[1:2],all.x=T)

q.wq$TP.int=with(q.wq,ave(TP,SITE,FUN=function(x) dat.interp(x)))
q.wq$TP.load=with(q.wq,Load.Calc.kg(Data.Value,TP.int))
q.wq$TN.int=with(q.wq,ave(TN,SITE,FUN=function(x) dat.interp(x)))
q.wq$TN.load=with(q.wq,Load.Calc.kg(Data.Value,TN.int))

##
q.wq.WY=ddply(q.wq,c("SITE","WY"),summarise,TFlow.acft=sum(cfs.to.acftd(Data.Value),na.rm=T),TPLoad=sum(TP.load,na.rm=T),TNLoad=sum(TN.load,na.rm=T))
q.wq.TPload.WY=cast(q.wq,WY~SITE,value="TP.load",fun.aggregate = function(x) sum(x,na.rm=T))
q.wq.TPload.WY$C43=with(q.wq.TPload.WY,ifelse(S79<S77,0,S79-S77))            
q.wq.TPload.WY=merge(q.wq.TPload.WY,rf.dat.WY[,c("WY","TRF.cm")],"WY")
q.wq.TPload.WY=merge(q.wq.TPload.WY,q.cre.dat.xtab.WY,"WY")
q.wq.TPload.WY=merge(q.wq.TPload.WY,lakeO.stg.WY,"WY")
q.wq.TPload.WY$S79.TP.roll=c(rep(NA,2),zoo::rollmean(q.wq.TPload.WY$S79,3))
q.wq.TPload.WY$S79.FWM=with(q.wq.TPload.WY,(S79/(Q.S79*4046.87260987425))*1e6)

# Data explore
x=boxplot(q.wq.TPload.WY$S79)
subset(q.wq.TPload.WY,S79%in%x$out)

plot(S79~WY,q.wq.TPload.WY,type="l",xlim=c(1980,2019))
with(subset(q.wq.TPload.WY,WY%in%cre.hurr$WY.start),points(WY,S79,pch=21,bg="red"))
par(new=T);plot(TRF.cm~WY,q.wq.TPload.WY,type="l",col="red",xlim=c(1980,2019))
par(new=T);plot(Q.C43~WY,q.wq.TPload.WY,type="l",col="blue",lty=2,xlim=c(1980,2019))

plot(S79~TRF.cm,q.wq.TPload.WY)
plot(S79~basin.ratio,q.wq.TPload.WY)
plot(S79.FWM~WY,q.wq.TPload.WY,type="l",xlim=c(1980,2019))
par(new=T);plot(Q.C43~WY,q.cre.dat.xtab.WY,type="l",col="red",xlim=c(1980,2019))

cols=wesanderson::wes_palette("Zissou1",length(q.wq.TPload.WY$WY),"continuous")
plot(S79~Q.C43,q.wq.TPload.WY)
with(q.wq.TPload.WY,points(Q.C43,S79,pch=21,bg=cols))

plot(S79~q1.stg,subset(q.wq.TPload.WY,!(WY%in%cre.hurr$WY.start)))
plot(S79~Q.C43,subset(q.wq.TPload.WY,!(WY%in%cre.hurr$WY.start)))

q.wq.TPload.WY.ex=subset(q.wq.TPload.WY,!(WY%in%cre.hurr$WY.start))
#mod.TP=lm(S79~Q.C43+Q.S77,subset(q.wq.TPload.WY,!(WY%in%c(1983,2006,2018))))
mod.TP=lm(sqrt(S79)~Q.C43+Q.S77+q1.stg,q.wq.TPload.WY.ex)
layout(matrix(1:4,2,2));plot(mod.TP)
dev.off()
summary(mod.TP)
shapiro.test(residuals(mod.TP))
acf(residuals(mod.TP))
gvlma::gvlma(mod.TP)
#AIC(mod.TP)

plot(S79~WY,q.wq.TPload.WY.ex,type="l")
mod.pred=predict(mod.TP,q.wq.TPload.WY.ex)
lines(q.wq.TPload.WY.ex$WY,mod.pred^2,col="red")


## Train vs test
# http://r-statistics.co/Linear-Regression.html
# using random 70% of data

set.seed(2)
tr.index=sample(1:nrow(q.wq.TPload.WY),nrow(q.wq.TPload.WY)*0.7)

q.wq.TPload.WY[tr.index,"index"]="train"
q.wq.TPload.WY[-tr.index,"index"]="test"
boxplot(S79~index,q.wq.TPload.WY)
kruskal.test(S79~index,q.wq.TPload.WY)

plot(S79~WY,q.wq.TPload.WY)
with(q.wq.TPload.WY[tr.index,],points(WY,S79,pch=21,bg="red"))

mod.TP=lm(S79~Q.C43,q.wq.TPload.WY[tr.index,])
layout(matrix(1:4,2,2));plot(mod.TP)
shapiro.test(residuals(mod.TP))
acf(mod.TP$residuals)
summary(mod.TP)
AIC(mod)
