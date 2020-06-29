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

library(mblm)
library(relaimpo)
library(car)

# GIS libraries 
library(rgdal)
library(rgeos)
library(raster)

## Paths
wd="C:/Julian_LaCie/_Github/LOSOM_WQ"

paths=paste0(wd,c("/Plots/WYModel/","/Export/","/Data/"))
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
ogrListLayers(paste0(GIS.path,"/SFER_GIS_Geodatabase.gdb"))

basins=spTransform(readOGR(paste0(GIS.path,"/AHED_release/AHED_20171102.gdb"),"WATERSHED"),utm17)
wmd.struct=spTransform(readOGR(paste0(GIS.path,"/AHED_release/AHED_20171102.gdb"),"STRUCTURE"),utm17)
shore=spTransform(readOGR(paste0(GIS.path,"/FWC"),"FWC_Shoreline"),utm17)
#shore=spTransform(readOGR(paste0(GIS.path,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Shoreline"),utm17)
canals=spTransform(readOGR(paste0(GIS.path,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Canals"),utm17)

library(tmap)
tmap_mode("view")
tm_shape(basins)+tm_polygons()

CRE.basins=c("WEST CALOOSAHATCHEE","EAST CALOOSAHATCHEE",
             "NICODEMUS SLOUGH NORTH","NICODEMUS SLOUGH SOUTH","HICPOCHEE NORTH",
             "S-4","LAKE OKEECHOBEE")#"CALOOSAHATCHEE ESTUARY"
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

cre.hurr=subset(cre.hurr,WY.start>=1981)
cre.hurr.trk=subset(hur.tck,Key%in%cre.hurr$Key)

plot(subset(basins,NAME%in%CRE.basins))
plot(subset(hur.tck,Key%in%cre.hurr$Key),add=T)

#cols=wesanderson::wes_palette("Zissou1",nrow(cre.hurr),"continuous")
#RColorBrewer::display.brewer.all(colorblindFriendly = T)
#cols=RColorBrewer::brewer.pal(nrow(cre.hurr),"PRGn")
#cols=viridis::viridis(nrow(cre.hurr))
cols=colorRampPalette(c("indianred1","forestgreen"))(nrow(cre.hurr))
bbox.lims=bbox(subset(basins,NAME%in%c(CRE.basins,"CALOOSAHATCHEE ESTUARY")))
#tiff(filename=paste0(plot.path,"CRE_HurricaneMap.tiff"),width=6.5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(1,1,1,1),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(1:2,1,2,byrow=F),widths=c(1,0.5))

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(basins,add=T,border="grey",lty=2,lwd=0.5,col=NA)
plot(subset(basins,NAME%in%CRE.basins),col=adjustcolor("grey",0.5),border="white",add=T)
plot(cre.buffer1,add=T,col=NA)
plot(canals,add=T,col="dodgerblue2",lwd=1.5)
plot(subset(wmd.struct,NAME%in%c("S79","S78","S77")),add=T,pch=21,bg="grey",cex=1.5)
plot(cre.hurr.trk[order(cre.hurr.trk$WY.start),],add=T,lwd=3,col=cols,lty=rep(c(1,2,3),3))
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4);box(lwd=1)

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA);
leg.text=c("C43 Basin","20 km Buffer","C43 Structures","Canals")
legend(0.5,0.85,legend=leg.text,
       pch=c(22,22,21,NA),
       lty=c(NA,NA,NA,1),lwd=c(1,1),
       col=c("white","black","black","dodgerblue2"),pt.bg=c(adjustcolor("grey",0.5),NA,"grey",NA),
       pt.cex=2,ncol=1,cex=0.8,bty="n",y.intersp=1.1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
leg.text=paste0(stringr::str_to_sentence(cre.hurr$Name)," (",cre.hurr$WY.start,")")
legend(0.5,0.3,legend=leg.text,
       pch=c(NA),
       lty=rep(c(1,2,3),3),lwd=3,
       col=cols,
       pt.cex=1,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title="Major Hurricane\nName (WY)",title.adj = 0)
dev.off()


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


# TP Model ----------------------------------------------------------------
q.wq.WY=ddply(q.wq,c("SITE","WY"),summarise,TFlow.acft=sum(cfs.to.acftd(Data.Value),na.rm=T),TPLoad=sum(TP.load,na.rm=T),TNLoad=sum(TN.load,na.rm=T))
q.wq.TPload.WY=cast(q.wq,WY~SITE,value="TP.load",fun.aggregate = function(x) sum(x,na.rm=T))
q.wq.TPload.WY$C43=with(q.wq.TPload.WY,ifelse(S79<S77,0,S79-S77))            
q.wq.TPload.WY=merge(q.wq.TPload.WY,rf.dat.WY[,c("WY","TRF.cm")],"WY")
q.wq.TPload.WY=merge(q.wq.TPload.WY,q.cre.dat.xtab.WY,"WY")
q.wq.TPload.WY=merge(q.wq.TPload.WY,lakeO.stg.WY,"WY")
q.wq.TPload.WY$S79.TP.roll=c(rep(NA,2),zoo::rollmean(q.wq.TPload.WY$S79,3))
q.wq.TPload.WY$S79.FWM=with(q.wq.TPload.WY,(S79/(Q.S79*1.233e6))*1e9)

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

plot(S79.FWM~TRF.cm,q.wq.TPload.WY)
plot(S79.FWM~basin.ratio,q.wq.TPload.WY)
plot(S79.FWM~Q.C43,q.wq.TPload.WY)

cols=wesanderson::wes_palette("Zissou1",length(q.wq.TPload.WY$WY),"continuous")
plot(S79~Q.C43,q.wq.TPload.WY)
with(q.wq.TPload.WY,points(Q.C43,S79,pch=21,bg=cols))

plot(S79~q1.stg,subset(q.wq.TPload.WY,!(WY%in%cre.hurr$WY.start)))
plot(S79~Q.C43,subset(q.wq.TPload.WY,!(WY%in%cre.hurr$WY.start)))

q.wq.TPload.WY.ex=subset(q.wq.TPload.WY,!(WY%in%cre.hurr$WY.start))
mod.TP=lm(sqrt(S79)~Q.C43+Q.S77+mean.stg,q.wq.TPload.WY.ex)
layout(matrix(1:4,2,2));plot(mod.TP)
dev.off()
summary(mod.TP)
shapiro.test(residuals(mod.TP))
acf(residuals(mod.TP))
gvlma::gvlma(mod.TP)

#AIC(mod.TP)

plot(S79~WY,q.wq.TPload.WY,type="l")
mod.pred=predict(mod.TP,q.wq.TPload.WY.ex)
points(q.wq.TPload.WY.ex$WY,mod.pred^2,col="red")
segments(q.wq.TPload.WY.ex$WY,mod.pred^2,q.wq.TPload.WY.ex$WY,q.wq.TPload.WY.ex$S79,col="red")

## Train vs test
# http://r-statistics.co/Linear-Regression.html
# using random 70% of data
set.seed(123)
tr.index=sample(1:nrow(q.wq.TPload.WY.ex),nrow(q.wq.TPload.WY.ex)*0.7)

q.wq.TPload.WY.ex[tr.index,"index"]="train"
q.wq.TPload.WY.ex[-tr.index,"index"]="test"
boxplot(S79~index,q.wq.TPload.WY.ex)
kruskal.test(S79~index,q.wq.TPload.WY.ex)

plot(S79~WY,q.wq.TPload.WY.ex)
with(q.wq.TPload.WY.ex[tr.index,],points(WY,S79,pch=21,bg="red"))

mod.TP=lm(sqrt(S79)~Q.C43+Q.S77+mean.stg,q.wq.TPload.WY.ex[tr.index,])
layout(matrix(1:4,2,2));plot(mod.TP)
shapiro.test(residuals(mod.TP))
acf(mod.TP$residuals)
summary(mod.TP)
car::vif(mod.TP)

gvlma::gvlma(mod.TP)
lmtest::bgtest(mod.TP)

boot <- boot.relimp(mod.TP, b = 1000, rank=T,diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE))
rslt=booteval.relimp(boot,sort=TRUE)
pick=which.max(rslt@level)
index <- sort(rslt@lmg, decreasing = T, index = T)$ix
xnames=rslt@namen[2:((length(rslt@namen) - 1) + 1)][index];xnames
xlabs=c(expression(paste("Q"["C43 Basin"])),"Mean Lake\nStage",expression(paste("Q"["S77"])))

ylim.val=c(0,1.05);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"CRE_WY_TPMod_relaimpo.tiff"),width=3,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,3,0.25,0.5),oma=c(0.5,1,0.25,0.5));
cols=wesanderson::wes_palette("Zissou1",length(xnames),"continuous")#nationalparkcolors::park_palette("Everglades",5)
x=barplot(t(rbind(rslt@lmg[index],rep(NA,length(xnames)))),col=rev(cols), ylim =ylim.val,axes=F)
text(x[1]+diff(x)/2,cumsum(rslt@lmg[index])-(cumsum(rslt@lmg[index])-c(0,cumsum(rslt@lmg[index])[1:2]))/2,xlabs,adj=0)
axis_fun(2,ymaj,ymin,ymaj*100);box(lwd=1)
mtext(side=2,line=2.5,"Percent of R\u00B2")
dev.off()

mod.pred=predict(mod.TP,q.wq.TPload.WY.ex[-tr.index,],interval="confidence")
actuals_preds <-data.frame(cbind(actuals=q.wq.TPload.WY.ex[-tr.index,"S79"],predicted=(mod.pred[,1])^2))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#min_max_accuracy
mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  

#tiff(filename=paste0(plot.path,"C43TPLoad_WY_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(1,6)*1e5;by.x=1e5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
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
axis_fun(1,xmaj,xmin,xmaj/1e5,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e5);box(lwd=1)
mtext(side=1,line=1.75,"Predicted TP Load (x10\u2074 kg WY\u207B\u00B9)")
mtext(side=2,line=2.25,"Actual TP Load (x10\u2074 kg WY\u207B\u00B9)")
legend("topleft",legend=c("Testing Dataset","1:1 line","Correlation"),
       pch=c(21,NA,NA),lty=c(NA,1,1),lwd=c(0.1,1,2),
       col=c(adjustcolor("black",0.5),"black","indianred1"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

kmodel.TP=data.frame()
kerrors.TP=data.frame()
k=10
set.seed(123)
for(i in 1:k){
  tr.index=sample(1:nrow(q.wq.TPload.WY.ex),nrow(q.wq.TPload.WY.ex)*0.7)
  
  train.dat=q.wq.TPload.WY.ex[tr.index,]
  test.dat=q.wq.TPload.WY.ex[-tr.index,]
  
  k.cv.mod=lm(sqrt(S79)~Q.C43+Q.S77+mean.stg,train.dat)
  summary(k.cv.mod)$r.squared
  
  k.cv.mod.pred=predict(k.cv.mod,test.dat)
  actuals_preds <-data.frame(cbind(actuals=as.numeric(test.dat$S79),predicted=(k.cv.mod.pred)^2))
  MMA.val=mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
  MAPE.val=mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
  
  actuals_preds$k=i
  kmodel.TP=rbind(kmodel.TP,actuals_preds)
  kerrors.TP=rbind(kerrors.TP,data.frame(k=i,R2.val=summary(k.cv.mod)$r.squared,RMSE=summary(k.cv.mod)$sigma,MMA=MMA.val,MAPE=MAPE.val))
  print(i)
}
apply(kerrors.TP,2,mean)

ylim.val=c(0,1);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.75,4.25)
plot(1:4,apply(kerrors.TP,2,mean)[2:5],ylim=c(0,1),type="n",axes=F,ylab=NA,xlab=NA,xlim=xlim.val)
arrows(1:4,apply(kerrors.TP[,2:5],2,min),1:4,apply(kerrors.TP[,2:5],2,max),angle=90,length=0.025,code=3,lwd=1)
points(1:4,apply(kerrors.TP,2,mean)[2:5],pch=21,bg="indianred1",lwd=0.1)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:4,1:4,c("R2","RMSE","MMA","MAPE"));box(lwd=1)
mtext(side=2,line=3,"Cross Validation Error")

plot(R2.val~k,kerrors.TP)
plot(MAPE~k,kerrors.TP)
plot(MMA~k,kerrors.TP)

cols=wesanderson::wes_palette("Zissou1",k,"continuous")
#tiff(filename=paste0(plot.path,"C43TPLoad_WY_kmodel.tiff"),width=4.5,height=2.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2.5,2.5,0.25,0.25));
xlim.val=c(0.5,6)*1e5;by.x=1.5e5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,kmodel.TP,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
for(i in 1:k){
  with(subset(kmodel.TP,k==i),points(predicted,actuals,pch=21,bg=adjustcolor(cols[i],0.5),col=cols[i],lwd=0.1,cex=0.75))
  mod.cor=mblm(actuals~predicted,subset(kmodel.TP,k==i))
  x.val=with(subset(kmodel.TP,k==i),seq(min(predicted),max(predicted),length.out = 25))
  mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
  lines(x.val,mod.cor.pred[,1],col=cols[i],lwd=1,lty=2)
}
axis_fun(1,xmaj,xmin,format(xmaj/1e5),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e5));box(lwd=1)
mtext(side=1,line=2.5,"Predicted TP Load\n(x10\u2074 kg WY\u207B\u00B9)")
mtext(side=2,line=2.25,"Actual TP Load\n(x10\u2074 kg WY\u207B\u00B9)")
dev.off()

library(caret)
kmodel.TP2=train(
  formula(mod.TP),
  q.wq.TPload.WY.ex,
  method="lm",
  trControl=trainControl(
    method="cv",
    number=10,
    verboseIter=T
  )
)
kmodel.TP2
print(kmodel.TP2)
kmodel.TP2$results
kmodel.TP2$finalModel
summary(kmodel.TP2)
shapiro.test(residuals(kmodel.TP2))

plot(S79~WY,q.wq.TPload.WY)
mod.TP.pred=predict(mod.TP,q.wq.TPload.WY,interval="confidence")
lines(q.wq.TPload.WY$WY,mod.TP.pred[,1]^2)

kfold.mod.pred=predict(kmodel.TP2,q.wq.TPload.WY,interval="confidence")
lines(q.wq.TPload.WY$WY,(kfold.mod.pred)^2,col="red",lwd=2)

# TN Model ----------------------------------------------------------------
q.wq.TNload.WY=cast(q.wq,WY~SITE,value="TN.load",fun.aggregate = function(x) sum(x,na.rm=T))
q.wq.TNload.WY$C43=with(q.wq.TNload.WY,ifelse(S79<S77,0,S79-S77))            
q.wq.TNload.WY=merge(q.wq.TNload.WY,rf.dat.WY[,c("WY","TRF.cm")],"WY")
q.wq.TNload.WY=merge(q.wq.TNload.WY,q.cre.dat.xtab.WY,"WY")
q.wq.TNload.WY=merge(q.wq.TNload.WY,lakeO.stg.WY,"WY")
q.wq.TNload.WY$S79.FWM=with(q.wq.TNload.WY,(S79/(Q.S79*1.233e6))*1e6)

q.wq.TNload.WY.ex=subset(q.wq.TNload.WY,!(WY%in%cre.hurr$WY.start))
#mod.TN=lm(S79~Q.C43+Q.S77,subset(q.wq.TNload.WY,!(WY%in%c(1983,2006,2018))))
mod.TN=lm(S79~Q.C43+Q.S77+mean.stg,q.wq.TNload.WY.ex)
layout(matrix(1:4,2,2));plot(mod.TN)
dev.off()
summary(mod.TN)
shapiro.test(residuals(mod.TN))
acf(residuals(mod.TN))
gvlma::gvlma(mod.TN)
lmtest::bgtest(mod.TN)

## Train vs test
set.seed(123)
tr.index=sample(1:nrow(q.wq.TNload.WY.ex),nrow(q.wq.TNload.WY.ex)*0.7)

q.wq.TNload.WY.ex[tr.index,"index"]="train"
q.wq.TNload.WY.ex[-tr.index,"index"]="test"
boxplot(S79~index,q.wq.TNload.WY.ex)
kruskal.test(S79~index,q.wq.TNload.WY.ex)

plot(S79~WY,q.wq.TNload.WY.ex)
with(q.wq.TNload.WY.ex[tr.index,],points(WY,S79,pch=21,bg="red"))

mod.TN=lm(S79~Q.C43+Q.S77+mean.stg,q.wq.TNload.WY.ex[tr.index,])
layout(matrix(1:4,2,2));plot(mod.TN)
shapiro.test(residuals(mod.TN))
acf(mod.TN$residuals)
summary(mod.TN)
car::vif(mod.TN)

gvlma::gvlma(mod.TN)
lmtest::bgtest(mod.TN)

boot <- boot.relimp(mod.TN, b = 1000, rank=T,diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE))
rslt=booteval.relimp(boot,sort=TRUE)
pick=which.max(rslt@level)
index <- sort(rslt@lmg, decreasing = T, index = T)$ix
xnames=rslt@namen[2:((length(rslt@namen) - 1) + 1)][index];xnames
xlabs=c(expression(paste("Q"["C43 Basin"])),expression(paste("Q"["S77"])),"Mean Lake\nStage")

ylim.val=c(0,1.05);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"CRE_WY_TNMod_relaimpo.tiff"),width=3,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,3,0.25,0.5),oma=c(0.5,1,0.25,0.5));
cols=wesanderson::wes_palette("Zissou1",length(xnames),"continuous")#nationalparkcolors::park_palette("Everglades",5)
x=barplot(t(rbind(rslt@lmg[index],rep(NA,length(xnames)))),col=rev(cols), ylim =ylim.val,axes=F)
text(x[1]+diff(x)/2,cumsum(rslt@lmg[index])-(cumsum(rslt@lmg[index])-c(0,cumsum(rslt@lmg[index])[1:2]))/2,xlabs,adj=0)
axis_fun(2,ymaj,ymin,ymaj*100);box(lwd=1)
mtext(side=2,line=2.5,"Percent of R\u00B2")
dev.off()

mod.pred=predict(mod.TN,q.wq.TNload.WY.ex[-tr.index,],interval="confidence")
actuals_preds <-data.frame(cbind(actuals=q.wq.TNload.WY.ex[-tr.index,"S79"],predicted=mod.pred[,1]))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#min_max_accuracy
mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) *100
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  *100

#png(filename=paste0(plot.path,"png/C43TNLoad_WY_ActualPred.png"),width=4,height=3,units="in",res=200,type="windows",bg="white")
#tiff(filename=paste0(plot.path,"C43TNLoad_WY_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2.5,2.5,0.25,0.25));
xlim.val=c(0.5,6)*1e6;by.x=2e6;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
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
axis_fun(1,xmaj,xmin,format(xmaj/1e6),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e6));box(lwd=1)
mtext(side=1,line=2.5,"Predicted TN Load\n(x10\u2075 kg WY\u207B\u00B9)")
mtext(side=2,line=2.25,"Actual TN Load\n(x10\u2075 kg WY\u207B\u00B9)")
legend("topleft",legend=c("Testing Dataset","1:1 line","Correlation"),
       pch=c(21,NA,NA),lty=c(NA,1,1),lwd=c(0.1,1,2),
       col=c(adjustcolor("black",0.5),"black","indianred1"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

kmodel.TN=data.frame()
kerrors.TN=data.frame()
k=10
set.seed(123)
for(i in 1:k){
  tr.index=sample(1:nrow(q.wq.TNload.WY.ex),nrow(q.wq.TNload.WY.ex)*0.7)
  
  train.dat=q.wq.TNload.WY.ex[tr.index,]
  test.dat=q.wq.TNload.WY.ex[-tr.index,]
  
  k.cv.mod=lm(S79~Q.C43+Q.S77+mean.stg,train.dat)
  summary(k.cv.mod)$r.squared
  
  k.cv.mod.pred=predict(k.cv.mod,test.dat)
  actuals_preds <-data.frame(cbind(actuals=as.numeric(test.dat$S79),predicted=k.cv.mod.pred))
  MMA.val=mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
  MAPE.val=mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
  
  actuals_preds$k=i
  kmodel.TN=rbind(kmodel.TN,actuals_preds)
  kerrors.TN=rbind(kerrors.TN,data.frame(k=i,R2.val=summary(k.cv.mod)$r.squared,RMSE=summary(k.cv.mod)$sigma,MMA=MMA.val,MAPE=MAPE.val))
  print(i)
}
apply(kerrors.TN,2,mean)

cols=wesanderson::wes_palette("Zissou1",k,"continuous")
#png(filename=paste0(plot.path,"png/C43TNLoad_WY_kmodel.png"),width=4,height=3,units="in",res=200,type="windows",bg="white")
#tiff(filename=paste0(plot.path,"C43TNLoad_WY_kmodel.tiff"),width=4.5,height=2.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2.5,2.5,0.25,0.25));
xlim.val=c(0.5,6)*1e6;by.x=2e6;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,kmodel.TN,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
for(i in 1:k){
  with(subset(kmodel.TN,k==i),points(predicted,actuals,pch=21,bg=adjustcolor(cols[i],0.5),col=cols[i],lwd=0.1,cex=0.75))
  mod.cor=mblm(actuals~predicted,subset(kmodel.TN,k==i))
  x.val=with(subset(kmodel.TN,k==i),seq(min(predicted),max(predicted),length.out = 25))
  mod.cor.pred=predict(mod.cor,data.frame(predicted=x.val),interval="confidence")
  lines(x.val,mod.cor.pred[,1],col=cols[i],lwd=1,lty=2)
}
axis_fun(1,xmaj,xmin,format(xmaj/1e6),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e6));box(lwd=1)
mtext(side=1,line=2.5,"Predicted TN Load\n(x10\u2075 kg WY\u207B\u00B9)")
mtext(side=2,line=2.25,"Actual TN Load\n(x10\u2075 kg WY\u207B\u00B9)")
dev.off()

## 

q.wq.TPload.WY$Pred.fit.TPLoad.kg=(predict(mod.TP,q.wq.TPload.WY,interval="confidence")[,1])^2
q.wq.TPload.WY$Pred.LCI.TPLoad.kg=(predict(mod.TP,q.wq.TPload.WY,interval="confidence")[,2])^2
q.wq.TPload.WY$Pred.UCI.TPLoad.kg=(predict(mod.TP,q.wq.TPload.WY,interval="confidence")[,3])^2

q.wq.TNload.WY$Pred.fit.TNLoad.kg=predict(mod.TN,q.wq.TPload.WY,interval="confidence")[,1]
q.wq.TNload.WY$Pred.LCI.TNLoad.kg=predict(mod.TN,q.wq.TPload.WY,interval="confidence")[,2]
q.wq.TNload.WY$Pred.UCI.TNLoad.kg=predict(mod.TN,q.wq.TPload.WY,interval="confidence")[,3]


#tiff(filename=paste0(plot.path,"C43_WY_ObsPredloads.tiff"),width=6.5,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,0.25,1),oma=c(2,1.5,0.75,0.25));
layout(matrix(1:2,1,2))
txt.cex=0.8
xlim.val=c(1,6)*1e5;by.x=1e5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.fit.TPLoad.kg~S79,q.wq.TPload.WY,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1.25,lty=2)
with(q.wq.TPload.WY,arrows(S79,Pred.LCI.TPLoad.kg,S79,Pred.UCI.TPLoad.kg,angle=90,code=3,length=0,col=adjustcolor("black",0.25)))
with(subset(q.wq.TPload.WY,WY%in%cre.hurr$WY.start),points(S79,Pred.fit.TPLoad.kg,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),cex=1.25,lwd=0.1))
with(subset(q.wq.TPload.WY,!(WY%in%cre.hurr$WY.start)),points(S79,Pred.fit.TPLoad.kg,pch=21,bg=adjustcolor("indianred1",0.5),col=adjustcolor("black",0.5),cex=1.25,lwd=0.1))
axis_fun(1,xmaj,xmin,xmaj/1e5,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e5);box(lwd=1)
mtext(side=1,line=1.75,"Predicted TP Load (x10\u2074 kg WY\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TP Load (x10\u2074 kg WY\u207B\u00B9)")
legend("topleft",legend=c("Hurricane Years"),pch=21,lty=NA,lwd=0.1,
       col="black",pt.bg=adjustcolor("dodgerblue1",0.5),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

xlim.val=c(0.5,6)*1e6;by.x=2e6;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(Pred.fit.TNLoad.kg~S79,q.wq.TNload.WY,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1.25,lty=2)
with(q.wq.TNload.WY,arrows(S79,Pred.LCI.TNLoad.kg,S79,Pred.UCI.TNLoad.kg,angle=90,code=3,length=0,col=adjustcolor("black",0.25)))
with(subset(q.wq.TNload.WY,WY%in%cre.hurr$WY.start),points(S79,Pred.fit.TNLoad.kg,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("black",0.5),cex=1.25,lwd=0.1))
with(subset(q.wq.TNload.WY,!(WY%in%cre.hurr$WY.start)),points(S79,Pred.fit.TNLoad.kg,pch=21,bg=adjustcolor("indianred1",0.5),col=adjustcolor("black",0.5),cex=1.25,lwd=0.1))
axis_fun(1,xmaj,xmin,format(xmaj/1e6),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj/1e6));box(lwd=1)
mtext(side=1,line=1.75,"Predicted TN Load (x10\u2075 kg WY\u207B\u00B9)")
mtext(side=2,line=2.5,"Actual TN Load (x10\u2075 kg WY\u207B\u00B9)")
dev.off()

#tiff(filename=paste0(plot.path,"C43_WYObsPredloads_annual.tiff"),width=4,height=5.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,0.25,1),oma=c(1,1.5,0.75,0.25));
layout(matrix(1:3,3,1,byrow=T),heights = c(1,1,0.2))
txt.cex=0.8

xlim.val=c(1982,2019);by.x=8;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,6.5)*1e5;by.y=2e5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(S79~WY,q.wq.TPload.WY,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
points(cre.hurr$WY.start,rep(ylim.val[2],nrow(cre.hurr)),pch=8,cex=0.8)
#abline(v=cre.hurr$WY.start,lty=2,col="grey50")
with(q.wq.TPload.WY,shaded.range(WY,Pred.LCI.TPLoad.kg,Pred.UCI.TPLoad.kg,"forestgreen",lty=0))
with(q.wq.TPload.WY,pt_line(WY,Pred.fit.TPLoad.kg,2,"forestgreen",1,21,"forestgreen"))
with(q.wq.TPload.WY,pt_line(WY,S79,2,"indianred1",1,21,"indianred1"))
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e5);box(lwd=1)
mtext(side=2,line=2.5,"TP Load (x10\u2074 kg WY\u207B\u00B9)")

ylim.val=c(0,6)*1e6;by.y=2e6;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(S79~WY,q.wq.TNload.WY,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
points(cre.hurr$WY.start,rep(ylim.val[2],nrow(cre.hurr)),pch=8,cex=0.8)
with(q.wq.TNload.WY,shaded.range(WY,Pred.LCI.TNLoad.kg,Pred.UCI.TNLoad.kg,"forestgreen",lty=0))
with(q.wq.TNload.WY,pt_line(WY,Pred.fit.TNLoad.kg,2,"forestgreen",1,21,"forestgreen"))
with(q.wq.TNload.WY,pt_line(WY,S79,2,"indianred1",1,21,"indianred1"))
abline(h=9086094*0.453592,lty=2,col="red");# Caloosahatchee TMDL 62-304.800(2), F.A.C.
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e6);box(lwd=1)
mtext(side=2,line=2.5,"TN Load (x10\u2075 kg WY\u207B\u00B9)")
mtext(side=1,line=1.75,"Florida Water Year (May - Apirl)",cex=txt.cex)

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA);
legend(0.5,-0.25,legend=c("Observed Load","Predicted Load \u00B1 95% CI","Hurricane Years"),
       pch=c(21,21,8),lty=c(NA,NA,NA),lwd=0.1,
       col=c("black","black","black"),pt.bg=c("indianred1","forestgreen",NA),
       pt.cex=1.5,ncol=3,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
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

RSM.lakeO.stg=ddply(RSM.lakeO.stg,c("WY"),summarise,
                        mean.stg=mean(STAGE,na.rm=T))

RSM.sites=c("S77","S78","S79")
RSM.q.cre.dat=data.frame()
for(i in 1:length(RSM.sites)){
  paths=paste0("/RSMBN/",RSM.sites[i],"/FLOW/01JAN1965 - 01JAN2016/1DAY/SIMULATED/")  
  tmp=data.frame(getFullTSC(dss_out,paths))
  tmp$Date=date.fun(rownames(tmp))
  tmp$SITE=RSM.sites[i]
  RSM.q.cre.dat=rbind(tmp,RSM.q.cre.dat)
  print(i)
}

RSM.q.cre.dat.xtab=cast(RSM.q.cre.dat,Date~SITE,value="FLOW",mean)
RSM.q.cre.dat.xtab$WY=WY(RSM.q.cre.dat.xtab$Date)
summary(RSM.q.cre.dat.xtab)

RSM.q.cre.dat.xtab.mon=ddply(RSM.q.cre.dat.xtab,c("WY"),summarise,Q.S77=sum(cfs.to.acftd(S77),na.rm=T),Q.S79=sum(cfs.to.acftd(S79),na.rm=T))
RSM.q.cre.dat.xtab.mon$Q.C43=with(RSM.q.cre.dat.xtab.mon,ifelse(Q.S79<Q.S77,0,Q.S79-Q.S77))

RSM.cre.hydro=merge(RSM.q.cre.dat.xtab.mon,RSM.lakeO.stg,c("WY"))
RSM.cre.hydro=subset(RSM.cre.hydro,WY%in%seq(1966,2016,1))

RSM.cre.hydro$TP.fit=(predict(mod.TP,RSM.cre.hydro,interval="confidence")[,1])^2
RSM.cre.hydro$TP.LCI=(predict(mod.TP,RSM.cre.hydro,interval="confidence")[,2])^2
RSM.cre.hydro$TP.UCI=(predict(mod.TP,RSM.cre.hydro,interval="confidence")[,3])^2
RSM.cre.hydro$S79.TPFWM=with(RSM.cre.hydro,(TP.fit/(Q.S79*1.233e6))*1e9)

RSM.cre.hydro$TN.fit=(predict(mod.TN,RSM.cre.hydro,interval="confidence")[,1])
RSM.cre.hydro$TN.LCI=(predict(mod.TN,RSM.cre.hydro,interval="confidence")[,2])
RSM.cre.hydro$TN.UCI=(predict(mod.TN,RSM.cre.hydro,interval="confidence")[,3])
RSM.cre.hydro$S79.TNFWM=with(RSM.cre.hydro,(TN.fit/(Q.S79*1.233e6))*1e6)


#tiff(filename=paste0(plot.path,"RSM_C43_WYPredloads_annual.tiff"),width=6.5,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,0.25,1),oma=c(2,2,0.75,0.25));
layout(matrix(1:4,2,2,byrow = T))
txt.cex=0.8
cols=adjustcolor("forestgreen",0.5)
xlim.val=c(1966,2016);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,6.5)*1e5;by.y=2e5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TP.fit~WY,RSM.cre.hydro,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.cre.hydro,pt_line(WY,TP.fit,1,cols,1,21,cols,pt.lwd=0.01))
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e5);box(lwd=1)
mtext(side=2,line=2.5,"TP Load (x10\u2074 kg WY\u207B\u00B9)")
legend("topleft",legend=c("LOSOM BN - LORS08"),
       pch=21,lty=NA,lwd=0.1,
       col="black",pt.bg=c("forestgreen"),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(100,200);by.y=25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(S79.TPFWM~WY,RSM.cre.hydro,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.cre.hydro,pt_line(WY,S79.TPFWM,1,cols,1,21,cols,pt.lwd=0.1))
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"TP FWM (\u03BCg L\u207B\u00B9)")

ylim.val=c(0,6)*1e6;by.y=2e6;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TN.fit~WY,RSM.cre.hydro,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.cre.hydro,pt_line(WY,TN.fit,1,cols,1,21,cols,pt.lwd=0.1))
abline(h=9086094*0.453592,lty=2,col="red");# Caloosahatchee TMDL 62-304.800(2), F.A.C.
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e6);box(lwd=1)
mtext(side=2,line=2.5,"TN Load (x10\u2075 kg WY\u207B\u00B9)")

ylim.val=c(1,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(S79.TNFWM~WY,RSM.cre.hydro,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA);
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(RSM.cre.hydro,pt_line(WY,S79.TNFWM,1,cols,1,21,cols,pt.lwd=0.1))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"TN FWM (mg L\u207B\u00B9)")
mtext(side=1,"Florida Water Year",outer=T,line=1)
dev.off()

vars=c("WY","Q.S77","Q.S79","Q.C43")
obs.q.sum=ddply(melt(q.cre.dat.xtab.WY[,vars],id.vars = "WY"),"variable",summarise,mean.Q=mean(value,na.rm=T),SE.Q=AnalystHelper::SE(value),N.Q=N.obs(value))
RSM.LORS.q.sum=ddply(melt(RSM.cre.hydro[,vars],id.vars = "WY"),"variable",summarise,mean.Q=mean(value,na.rm=T),SE.Q=AnalystHelper::SE(value),N.Q=N.obs(value))

tmp=rbind(t(obs.q.sum[match(obs.q.sum$variable,c("Q.S77","Q.C43","Q.S79")),"mean.Q"]),
      t(RSM.LORS.q.sum[match(RSM.LORS.q.sum$variable,c("Q.S77","Q.C43","Q.S79")),"mean.Q"]))
tmp.SE=rbind(t(obs.q.sum[match(obs.q.sum$variable,c("Q.S77","Q.C43","Q.S79")),"SE.Q"]),
          t(RSM.LORS.q.sum[match(RSM.LORS.q.sum$variable,c("Q.S77","Q.C43","Q.S79")),"SE.Q"]))
#
ylim.val=c(0,170e4);by.y=40e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylabs=c(expression(paste("Q"["S77"])),expression(paste("Q"["C43 Basin"])),expression(paste("Q"["S79"])))
#tiff(filename=paste0(plot.path,"RSM_Obs_annualQ.tiff"),width=6,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,2,0.25,0.25),oma=c(2,2,0.75,0.25));
layout(matrix(1:2,1,2,byrow = T),widths=c(1,0.5))

x=barplot(tmp,beside=T,horiz=T,xlim=ylim.val,axes=F,col=adjustcolor(c("indianred1","dodgerblue1"),0.5))
arrows(tmp-tmp.SE,x,tmp+tmp.SE,x,length=0.05,angle=90,code=3)
axis_fun(2,x[1,]+(diff(x)/2),x[1,]+(diff(x)/2),ylabs)
axis_fun(1,ymaj,ymin,ymaj/1e4,line=-0.5);box(lwd=1)
mtext(side=1,line=1.75,"Discharge (x1000 Acre-Ft WY\u207B\u00B9)")

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA);
legend(0.30,0.5,legend=c("RSM-BN LORS08\n(WY1966 - 2016)","Observed\n(WY1980 - 2019)"),
       pch=22,lty=NA,lwd=0.1,
       col="black",pt.bg=adjustcolor(c("dodgerblue1","indianred1"),0.5),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1.25,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title="Scenario\n(Mean \u00B1 SE)",title.adj = 0)
dev.off()

plot(Q.S77~WY,RSM.cre.hydro,xlim=c(1966,2019),type="l")
with(q.cre.dat.xtab.WY,lines(WY,Q.S77,col="red"))
