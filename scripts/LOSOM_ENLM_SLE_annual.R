## 
## LOSOM Water Quality Subteam
## Estuary Nutrient Load Model
## St Lucie Estuary
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

#
nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##
dates=as.Date(c("1978-05-01","2019-04-30"))

# Discharge  --------------------------------------------------------------
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
q.sle.dat$WY=WY(q.sle.dat$Date.EST)

# q.sle.dat$delta.day=with(q.sle.dat,c(NA,diff(Date.EST)))
# range(subset(q.sle.dat,SITE=="S308")$delta.day,na.rm = T)
# range(subset(q.sle.dat,SITE=="S80")$delta.day,na.rm = T)
# subset(q.sle.dat,SITE=="S308"&delta.day>1)
# subset(q.sle.dat,SITE=="S80"&delta.day!=1)

q.sle.dat.xtab=cast(q.sle.dat,Date.EST~SITE,value="Data.Value",mean)
q.sle.dat.xtab$month=format(q.sle.dat.xtab$Date,"%m")
q.sle.dat.xtab$CY=format(q.sle.dat.xtab$Date,"%Y")
q.sle.dat.xtab$monCY=with(q.sle.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
q.sle.dat.xtab$WY=WY(q.sle.dat.xtab$Date.EST)

#Positive discharge only
q.sle.dat.xtab$S308=with(q.sle.dat.xtab,ifelse(S308<0,0,S308))
q.sle.dat.xtab$S80=with(q.sle.dat.xtab,ifelse(S80<0,0,S80))
q.sle.dat.xtab$C44=with(q.sle.dat.xtab,ifelse(S80<S308,0,S80-S308))

q.sle.dat.xtab.WY=ddply(q.sle.dat.xtab,c("WY"),summarise,Q.S308=sum(cfs.to.acftd(S308),na.rm=T),Q.S80=sum(cfs.to.acftd(S80),na.rm=T))
q.sle.dat.xtab.WY$Q.C44=with(q.sle.dat.xtab.WY,ifelse(Q.S80<Q.S308,0,Q.S80-Q.S308))
q.sle.dat.xtab.WY$basin.ratio=with(q.sle.dat.xtab.WY,Q.C44/Q.S80)

# CRE WQ ------------------------------------------------------------------
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7),param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp"))
wq.sites=c("C44S80","S308C")
wq.dat=DBHYDRO_WQ(dates[1],dates[2],wq.sites,params$Test.Number)
wq.dat=merge(wq.dat,params,"Test.Number")
unique(wq.dat$Collection.Method)
wq.dat=subset(wq.dat,Collection.Method=="G")

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
C44Basin_rainfall.sites=data.frame(DBKEY=c(c("06239","06119","16588"),c("05849","K8637"),c("06075","06237","16618"),"06122",c("SX445","VM862")),
                                   SITE=c(rep("S308R",3),rep("S135R",2),rep("S80R",3),"PRATT",rep("ACRA2",2)))
C44Basin_rainfall.sites$BASIN="C44"

rf.dat=data.frame()
pb=txtProgressBar(1,nrow(C44Basin_rainfall.sites),1,style=3)
for(i in 1:nrow(C44Basin_rainfall.sites)){
  tmp=DBHYDRO_daily(dates[1],dates[2],C44Basin_rainfall.sites$DBKEY[i])
  tmp$DBKEY=as.character(C44Basin_rainfall.sites$DBKEY[i])
  rf.dat=rbind(tmp,rf.dat)
  setTxtProgressBar(pb, i)
}
rf.dat=merge(rf.dat,C44Basin_rainfall.sites,"DBKEY")
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

# library(tmap)
# tmap_mode("view")
# tm_shape(basins)+tm_polygons()
# 
SLE.basins=c("C-44","LAKE OKEECHOBEE")#UPPER ST. LUCIE ESTUARY
# plot(subset(basins,NAME%in%SLE.basins))

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

#SLE boundary and buffer
plot(subset(basins,NAME%in%SLE.basins))
sle.buffer1=gBuffer(subset(basins,NAME%in%SLE.basins),width=20*1000);
plot(sle.buffer1,add=T)
sle.buffer1=SpatialPolygonsDataFrame(sle.buffer1,data.frame(row.names = "buffer",width.km=20,area="SLE"))
sle.buffer1=spTransform(sle.buffer1,utm17)

sle.hurr=over(sle.buffer1,hur.tck,returnList = T,byid=T)[[1]];#select only hurricanes that cross the 200 km buffer. 
sle.hurr=sle.hurr[order(sle.hurr$WY.start,sle.hurr$Key),]

sle.hurr=subset(sle.hurr,WY.start>=1980)
sle.hurr.trk=subset(hur.tck,Key%in%sle.hurr$Key)

plot(subset(basins,NAME%in%SLE.basins))
plot(subset(hur.tck,Key%in%sle.hurr$Key),add=T)

#cols=wesanderson::wes_palette("Zissou1",nrow(cre.hurr),"continuous")
#RColorBrewer::display.brewer.all(colorblindFriendly = T)
#cols=RColorBrewer::brewer.pal(nrow(cre.hurr),"PRGn")
#cols=viridis::viridis(nrow(cre.hurr))
cols=colorRampPalette(c("indianred1","forestgreen"))(nrow(sle.hurr))
#bbox.lims=bbox(subset(basins,NAME%in%c(SLE.basins,"UPPER ST. LUCIE ESTUARY")))
bbox.lims=bbox(sle.buffer1)
#tiff(filename=paste0(plot.path,"SLE_HurricaneMap.tiff"),width=6.5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(1,1,1,1),mar=c(0.1,0.1,0.1,0.1),xpd=F)
layout(matrix(1:2,1,2,byrow=F),widths=c(1,0.5))

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(basins,add=T,border="grey",lty=2,lwd=0.5,col=NA)
plot(subset(basins,NAME%in%SLE.basins),col=adjustcolor("grey",0.5),border="white",add=T)
plot(sle.buffer1,add=T,col=NA)
plot(canals,add=T,col="dodgerblue2",lwd=1.5)
plot(subset(wmd.struct,NAME%in%c("S308","S80")),add=T,pch=21,bg="grey",cex=1.5)
plot(sle.hurr.trk[order(sle.hurr.trk$WY.start),],add=T,lwd=3,col=cols,lty=rep(c(1,2,3),nrow(sle.hurr)))
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4);box(lwd=1)

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA);
leg.text=c("C44 Basin","20 km Buffer","C44 Structures","Canals")
legend(0.5,0.85,legend=leg.text,
       pch=c(22,22,21,NA),
       lty=c(NA,NA,NA,1),lwd=c(1,1),
       col=c("white","black","black","dodgerblue2"),pt.bg=c(adjustcolor("grey",0.5),NA,"grey",NA),
       pt.cex=2,ncol=1,cex=0.8,bty="n",y.intersp=1.1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
leg.text=paste0(stringr::str_to_sentence(sle.hurr$Name)," (",sle.hurr$WY.start,")")
legend(0.5,0.3,legend=leg.text,
       pch=c(NA),
       lty=rep(c(1,2,3),3),lwd=3,
       col=cols,
       pt.cex=1,ncol=1,cex=0.7,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title="Major Hurricane\nName (WY)",title.adj = 0)
dev.off()

# write shp files
# writeOGR(basins,paste0(data.path,"GIS"),"basins",driver="ESRI Shapefile")
# writeOGR(canals,paste0(data.path,"GIS"),"canals",driver="ESRI Shapefile")
# writeOGR(wmd.struct,paste0(data.path,"GIS"),"wmd_struct",driver="ESRI Shapefile")
# writeOGR(shore,paste0(data.path,"GIS"),"FL_Shore",driver="ESRI Shapefile")
# writeOGR(hur.tck,paste0(data.path,"GIS"),"HURDAT_AL",driver="ESRI Shapefile")

# Load --------------------------------------------------------------------
head(q.sle.dat)
unique(q.sle.dat$SITE)
range(q.sle.dat$Data.Value,na.rm=T)
q.sle.dat$Data.Value=with(q.sle.dat,ifelse(Data.Value<0,NA,Data.Value))
head(wq.dat.xtab)
unique(wq.dat.xtab$Station.ID)
head(subset(wq.dat.xtab,Station.ID=="S308C"))
head(subset(wq.dat.xtab,Station.ID=="C44S80"))

wq.dat.xtab=merge(wq.dat.xtab,data.frame(Station.ID=wq.sites,SITE=c("S80","S308")),"Station.ID")
q.vars=c("SITE","Date.EST","WY","Data.Value")
wq.vars=c("SITE","Date.EST","TP","TN")
q.wq=merge(subset(q.sle.dat[,q.vars],WY>1981),
           wq.dat.xtab[,wq.vars],wq.vars[1:2],all.x=T)

#q.wq$TP=with(q.wq, ifelse(Data.Value==0|is.na(Data.Value)==T,NA,TP))
q.wq$TP.int=with(q.wq,ave(TP,SITE,FUN=function(x) dat.interp(x)))
q.wq$TP.load=with(q.wq,Load.Calc.kg(Data.Value,TP.int))
#q.wq$TN=with(q.wq, ifelse(Data.Value==0|is.na(Data.Value)==T,NA,TN))
q.wq$TN.int=with(q.wq,ave(TN,SITE,FUN=function(x) dat.interp(x)))
q.wq$TN.load=with(q.wq,Load.Calc.kg(Data.Value,TN.int))

# TP Model ----------------------------------------------------------------
q.wq.WY=ddply(q.wq,c("SITE","WY"),summarise,TFlow.acft=sum(cfs.to.acftd(Data.Value),na.rm=T),TPLoad=sum(TP.load,na.rm=T),TNLoad=sum(TN.load,na.rm=T))
q.wq.TPload.WY=cast(q.wq,WY~SITE,value="TP.load",fun.aggregate = function(x) sum(x,na.rm=T))
q.wq.TPload.WY$C44=with(q.wq.TPload.WY,ifelse(S80<S308,0,S80-S308))            
q.wq.TPload.WY=merge(q.wq.TPload.WY,q.sle.dat.xtab.WY,"WY")
q.wq.TPload.WY=merge(q.wq.TPload.WY,lakeO.stg.WY,"WY")

q.wq.TPload.WY$S80.FWM=with(q.wq.TPload.WY,(S80/(Q.S80*1.233e6))*1e9)
q.wq.TPload.WY=subset(q.wq.TPload.WY,Q.S80>120) # excluded extreme drought/no flow years 2008 and 2012

# Data explore
x=boxplot(q.wq.TPload.WY$S80)
subset(q.wq.TPload.WY,S80%in%x$out)

plot(S80~WY,q.wq.TPload.WY,type="l",xlim=c(1980,2019))
with(subset(q.wq.TPload.WY,WY%in%sle.hurr$WY.start),points(WY,S80,pch=21,bg="red"))
par(new=T);plot(Q.C44~WY,q.wq.TPload.WY,type="l",col="blue",lty=2,xlim=c(1980,2019))

plot(S80~basin.ratio,q.wq.TPload.WY)
plot(S80.FWM~WY,q.wq.TPload.WY,type="l",xlim=c(1980,2019))
par(new=T);plot(Q.C44~WY,q.sle.dat.xtab.WY,type="l",col="red",xlim=c(1980,2019))

q.wq.TPload.WY.ex=subset(q.wq.TPload.WY,!(WY%in%unique(sle.hurr$WY.start)))
mod.TP=lm(S80~Q.C44+Q.S308+mean.stg,q.wq.TPload.WY.ex)
layout(matrix(1:4,2,2));plot(mod.TP)
dev.off()
summary(mod.TP)
shapiro.test(residuals(mod.TP))
acf(residuals(mod.TP))
gvlma::gvlma(mod.TP)

plot(S80~WY,q.wq.TPload.WY,type="l")
mod.pred=predict(mod.TP,q.wq.TPload.WY.ex)
points(q.wq.TPload.WY.ex$WY,mod.pred,col="red")
segments(q.wq.TPload.WY.ex$WY,mod.pred,q.wq.TPload.WY.ex$WY,q.wq.TPload.WY.ex$S80,col="red")

## Train vs test
# http://r-statistics.co/Linear-Regression.html
# using random 70% of data
set.seed(1234)
tr.index=sample(1:nrow(q.wq.TPload.WY.ex),nrow(q.wq.TPload.WY.ex)*0.7)

q.wq.TPload.WY.ex[tr.index,"index"]="train"
q.wq.TPload.WY.ex[-tr.index,"index"]="test"
boxplot(S80~index,q.wq.TPload.WY.ex)
kruskal.test(S80~index,q.wq.TPload.WY.ex)

plot(S80~WY,q.wq.TPload.WY.ex)
with(q.wq.TPload.WY.ex[tr.index,],points(WY,S80,pch=21,bg="red"))

mod.TP=lm(sqrt(S80)~Q.C44+Q.S308+mean.stg,q.wq.TPload.WY.ex[tr.index,])
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
xlabs=c(expression(paste("Q"["S308"])),"Mean Lake\nStage",expression(paste("Q"["C44 Basin"])))

ylim.val=c(0,1.05);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"SLE_WY_TPMod_relaimpo.tiff"),width=3,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,3,0.25,0.5),oma=c(0.5,1,0.25,0.5));
cols=wesanderson::wes_palette("Zissou1",length(xnames),"continuous")#nationalparkcolors::park_palette("Everglades",5)
x=barplot(t(rbind(rslt@lmg[index],rep(NA,length(xnames)))),col=rev(cols), ylim =ylim.val,axes=F)
text(x[1]+diff(x)/2,cumsum(rslt@lmg[index])-(cumsum(rslt@lmg[index])-c(0,cumsum(rslt@lmg[index])[1:2]))/2,xlabs,adj=0)
axis_fun(2,ymaj,ymin,ymaj*100);box(lwd=1)
mtext(side=2,line=2.5,"Percent of R\u00B2")
dev.off()

mod.pred=predict(mod.TP,q.wq.TPload.WY.ex[-tr.index,],interval="confidence")^2
actuals_preds <-data.frame(cbind(actuals=q.wq.TPload.WY.ex[-tr.index,"S80"],predicted=mod.pred[,1]))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#min_max_accuracy
mean(apply(actuals_preds, 1, min,na.rm=T) / apply(actuals_preds, 1, max,na.rm=T))*100
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals,na.rm=T)*100

#tiff(filename=paste0(plot.path,"C44TPLoad_WY_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(0,3.1)*1e5;by.x=1e5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
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
set.seed(1234)
for(i in 1:k){
  tr.index=sample(1:nrow(q.wq.TPload.WY.ex),nrow(q.wq.TPload.WY.ex)*0.7)
  
  train.dat=q.wq.TPload.WY.ex[tr.index,]
  test.dat=q.wq.TPload.WY.ex[-tr.index,]
  
  k.cv.mod=lm(sqrt(S80)~Q.C44+Q.S308+mean.stg,train.dat)
  summary(k.cv.mod)$r.squared
  
  k.cv.mod.pred=predict(k.cv.mod,test.dat)^2
  actuals_preds <-data.frame(cbind(actuals=as.numeric(test.dat$S80),predicted=k.cv.mod.pred))
  MMA.val=mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
  MAPE.val=mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
  
  actuals_preds$k=i
  kmodel.TP=rbind(kmodel.TP,actuals_preds)
  kerrors.TP=rbind(kerrors.TP,data.frame(k=i,R2.val=summary(k.cv.mod)$r.squared,RMSE=summary(k.cv.mod)$sigma,MMA=MMA.val,MAPE=MAPE.val))
  print(i)
}
apply(kerrors.TP,2,mean,na.rm=T)

ylim.val=c(0,1);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0.75,4.25)
plot(1:4,apply(kerrors.TP,2,mean)[2:5],ylim=c(0,1),type="n",axes=F,ylab=NA,xlab=NA,xlim=xlim.val)
arrows(1:4,apply(kerrors.TP[,2:5],2,min),1:4,apply(kerrors.TP[,2:5],2,max),angle=90,length=0.025,code=3,lwd=1)
points(1:4,apply(kerrors.TP,2,mean)[2:5],pch=21,bg="indianred1",lwd=0.1)
axis_fun(2,ymaj,ymin,format(ymaj))
axis_fun(1,1:4,1:4,c("R2","RMSE","MMA","MAPE"));box(lwd=1)
mtext(side=2,line=3,"Cross Validation Error")


cols=wesanderson::wes_palette("Zissou1",k,"continuous")
#tiff(filename=paste0(plot.path,"C44TPLoad_WY_kmodel.tiff"),width=4.5,height=2.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2.5,2.5,0.25,0.25));
xlim.val=c(0,4)*1e5;by.x=1e5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=xlim.val;by.y=by.x;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(actuals~predicted,kmodel.TP,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(0,1,lwd=1)
for(i in 1:k){
  with(subset(kmodel.TP,k==i),points(predicted,actuals,pch=21,bg=adjustcolor(cols[i],0.5),col=cols[i],lwd=0.1,cex=0.75))
  mod.cor=mblm(actuals~predicted,subset(kmodel.TP,k==i))
  x.val=with(subset(kmodel.TP,k==i),seq(min(predicted,na.rm=T),max(predicted,na.rm=T),length.out = 25))
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

plot(S80~WY,q.wq.TPload.WY)
mod.TP.pred=predict(mod.TP,q.wq.TPload.WY,interval="confidence")^2
lines(q.wq.TPload.WY$WY,mod.TP.pred[,1])
kfold.mod.pred=predict(kmodel.TP2,q.wq.TPload.WY,interval="confidence")^2
lines(q.wq.TPload.WY$WY,(kfold.mod.pred),col="red",lwd=2)

# TN Model ----------------------------------------------------------------
q.wq.TNload.WY=cast(q.wq,WY~SITE,value="TN.load",fun.aggregate = function(x) sum(x,na.rm=T))
q.wq.TNload.WY$C44=with(q.wq.TNload.WY,ifelse(S80<S308,0,S80-S308))            
q.wq.TNload.WY=merge(q.wq.TNload.WY,q.sle.dat.xtab.WY,"WY")
q.wq.TNload.WY=merge(q.wq.TNload.WY,lakeO.stg.WY,"WY")
q.wq.TNload.WY$S80.FWM=with(q.wq.TNload.WY,(S80/(Q.S80*1.233e6))*1e6)
q.wq.TNload.WY=subset(q.wq.TNload.WY,Q.S80>120) # excluded extreme drought/no flow years 2008 and 2012

q.wq.TNload.WY.ex=subset(q.wq.TNload.WY,!(WY%in%sle.hurr$WY.start))
mod.TN=lm(S80~Q.S308+basin.ratio,q.wq.TNload.WY.ex)
layout(matrix(1:4,2,2));plot(mod.TN)
dev.off()
summary(mod.TN)
shapiro.test(residuals(mod.TN));hist(residuals(mod.TN))
acf(residuals(mod.TN))
car::vif(mod.TN)
gvlma::gvlma(mod.TN)
lmtest::bgtest(mod.TN)


## Train vs test
set.seed(1234)
tr.index=sample(1:nrow(q.wq.TNload.WY.ex),nrow(q.wq.TNload.WY.ex)*0.7)

q.wq.TNload.WY.ex[tr.index,"index"]="train"
q.wq.TNload.WY.ex[-tr.index,"index"]="test"
boxplot(S80~index,q.wq.TNload.WY.ex)
kruskal.test(S80~index,q.wq.TNload.WY.ex)

plot(S80~WY,q.wq.TNload.WY.ex)
with(q.wq.TNload.WY.ex[tr.index,],points(WY,S80,pch=21,bg="red"))

mod.TN=lm(sqrt(S80)~Q.S308+Q.C44+q1.stg,q.wq.TNload.WY.ex[tr.index,])
layout(matrix(1:4,2,2));plot(mod.TN)
dev.off()
summary(mod.TN)
shapiro.test(residuals(mod.TN));hist(residuals(mod.TN))
acf(residuals(mod.TN))
car::vif(mod.TN)
gvlma::gvlma(mod.TN)
lmtest::bgtest(mod.TN)

boot <- boot.relimp(mod.TN, b = 1000, rank=T,diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE))
rslt=booteval.relimp(boot)
pick=which.max(rslt@level)
index <- sort(rslt@lmg, decreasing = T, index = T)$ix
xnames=rslt@namen[2:((length(rslt@namen) - 1) + 1)][index];xnames
xlabs=c(expression(paste("Q"["S308"])),expression(paste("Q"["C44 Basin"])),expression(paste("q"["1"],"Stage")))

ylim.val=c(0,1.05);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"SLE_WY_TNMod_relaimpo.tiff"),width=3,height=3.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,3,0.25,0.5),oma=c(0.5,1,0.25,0.5));
cols=wesanderson::wes_palette("Zissou1",length(xnames),"continuous")#nationalparkcolors::park_palette("Everglades",5)
x=barplot(t(rbind(rslt@lmg[index],rep(NA,length(xnames)))),col=rev(cols), ylim =ylim.val,axes=F)
text(x[1]+diff(x)/2,cumsum(rslt@lmg[index])-(cumsum(rslt@lmg[index])-c(0,cumsum(rslt@lmg[index])[1:2]))/2,xlabs,adj=0)
axis_fun(2,ymaj,ymin,ymaj*100);box(lwd=1)
mtext(side=2,line=2.5,"Percent of R\u00B2")
dev.off()

mod.pred=predict(mod.TN,q.wq.TNload.WY.ex[-tr.index,],interval="confidence")^2
actuals_preds <-data.frame(cbind(actuals=q.wq.TNload.WY.ex[-tr.index,"S80"],predicted=mod.pred[,1]))
plot(actuals~predicted,actuals_preds);abline(0,1)
cor.test(actuals_preds$actuals,actuals_preds$predicted,method="spearman")

#min_max_accuracy
mean(apply(actuals_preds, 1, min,na.rm=T) / apply(actuals_preds, 1, max,na.rm=T))
#Mean Absolute Percent Error
mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals,na.rm=T)

#tiff(filename=paste0(plot.path,"C44TNLoad_WY_ActualPred.tiff"),width=4,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.75,0.25,0.25),oma=c(2,2,0.75,0.5));
xlim.val=c(0,3.0)*1e6;by.x=1e6;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
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
mtext(side=1,line=1.75,"Predicted TN Load (x10\u2075 kg WY\u207B\u00B9)")
mtext(side=2,line=2.25,"Actual TN Load (x10\u2075 kg WY\u207B\u00B9)")
legend("topleft",legend=c("Testing Dataset","1:1 line","Correlation"),
       pch=c(21,NA,NA),lty=c(NA,1,1),lwd=c(0.1,1,2),
       col=c(adjustcolor("black",0.5),"black","indianred1"),
       pt.bg=c(adjustcolor("dodgerblue1",0.5),NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

kmodel.TN=data.frame()
kerrors.TN=data.frame()
k=10
set.seed(1234)
for(i in 1:k){
  tr.index=sample(1:nrow(q.wq.TNload.WY.ex),nrow(q.wq.TNload.WY.ex)*0.7)
  
  train.dat=q.wq.TNload.WY.ex[tr.index,]
  test.dat=q.wq.TNload.WY.ex[-tr.index,]
  
  k.cv.mod=lm(sqrt(S80)~Q.S308+Q.C44+q1.stg,train.dat)
  summary(k.cv.mod)$r.squared
  
  k.cv.mod.pred=predict(k.cv.mod,test.dat)^2
  actuals_preds <-data.frame(cbind(actuals=as.numeric(test.dat$S80),predicted=k.cv.mod.pred))
  MMA.val=mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max)) 
  MAPE.val=mean(abs((actuals_preds$predicted - actuals_preds$actuals))/actuals_preds$actuals)  
  
  actuals_preds$k=i
  kmodel.TN=rbind(kmodel.TN,actuals_preds)
  kerrors.TN=rbind(kerrors.TN,data.frame(k=i,R2.val=summary(k.cv.mod)$r.squared,RMSE=summary(k.cv.mod)$sigma,MMA=MMA.val,MAPE=MAPE.val))
  print(i)
}
apply(kerrors.TN,2,mean)