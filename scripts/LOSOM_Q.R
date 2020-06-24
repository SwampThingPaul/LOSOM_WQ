## 
## LOSOM 
## Savings Clause Evaluation
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

getOption('timeout')
options(timeout=60*20)

#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)

library(openxlsx)

#Paths
wd="C:/Julian_LaCie/_Github/LOSOM_WQ"

paths=paste0(wd,c("/Plots/LakeO/","/Export/","/Data/"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

##
dates=as.Date(c("2009-05-01","2020-04-30"))

# Lake Okeechobee Discharge -----------------------------------------------
data.path2="C:/Julian_LaCie/Work/LakeOkeechobee/Data/"
flow.dbkeys=read.xlsx(paste0(data.path2,"LakeO_DBKEYS_V3.xlsx"),sheet=1)
# Double check dbhydro for dbkey start and end date
# paste(flow.dbkeys$DBKEY,collapse="/")
# flow.dbkeys=subset(flow.dbkeys,!(DBKEY%in%c("P0794","P1020","90880")));

# LO.q=data.frame()
# pb=txtProgressBar(1,nrow(flow.dbkeys),1,style=3)
# for(i in 1:nrow(flow.dbkeys)){
#     tmp=DBHYDRO_daily(dates[1],dates[2],flow.dbkeys$DBKEY[i])
#     tmp$DBKEY=as.character(flow.dbkeys$DBKEY[i])
#     LO.q=rbind(tmp,LO.q)
#     setTxtProgressBar(pb, i)
# }

LO.q=read.csv(paste0(data.path,"20200602_LakeOQ.csv"),skip=64)
LO.q$DBKEY=as.character(LO.q$DBKEY)
LO.q=merge(LO.q,flow.dbkeys,"DBKEY")
LO.q$Date=with(LO.q,as.POSIXct(as.character(Daily.Date),format="%d-%b-%Y",tz="America/New_York"))
LO.q$Date.EST=date.fun(LO.q$Date)
LO.q$WY=WY(LO.q$Date.EST)

flow.xtab=data.frame(cast(LO.q,Date.EST+WY+STRUCT+ALIAS+Inflow_direction+Outflow+Basin+WQSite~Priority,value="Data.Value",fun.aggregate=function(x) ifelse(sum(x,na.rm=T)==0,NA,sum(x,na.rm=T))))
flow.xtab$fflow.cfs=with(flow.xtab,ifelse(is.na(P1),P2,P1));#if only two priorities

flow.xtab$fflow.cfs=with(flow.xtab,fflow.cfs*Inflow_direction)#all inflows are positive and all negative values are outflow
flow.xtab$direct=with(flow.xtab,ifelse(fflow.cfs<0,"Outflow","Inflow"))
flow.xtab$month=as.numeric(format(flow.xtab$Date,"%m"))
flow.xtab$CY=as.numeric(format(flow.xtab$Date,"%Y"))

mon.seq=data.frame(Date.EST=date.fun(seq(dates[1],dates[2],"1 months")))
mon.seq$month=as.numeric(format(mon.seq$Date.EST,"%m"))
mon.seq$CY=as.numeric(format(mon.seq$Date.EST,"%Y"))

flow.mon.sum=cast(flow.xtab,STRUCT+ALIAS+Basin+WQSite+WY+CY+month~direct,value="fflow.cfs",fun.aggregate=function(x) sum(abs(x),na.rm=T))
flow.mon.sum=flow.mon.sum[,c("STRUCT", "ALIAS", "Basin", "WQSite", "WY", "CY", "month","Inflow", "Outflow")]
##
##
# Estimate flow for L31E, HP7 & Inflow 1,2,3
L61E_HP7=data.frame()
for(i in 1:nrow(mon.seq)){
  tmp.dat=subset(flow.mon.sum,month==mon.seq$month[i]&CY==mon.seq$CY[i])
  C41H78=subset(tmp.dat,ALIAS=="C41H78_I")
  G76=subset(tmp.dat,ALIAS=="G76_C")
  G207=subset(tmp.dat,ALIAS=="G207")
  S71=subset(tmp.dat,ALIAS=="S71_S")
  # if(nrow(C41H78)==0){
  #   L61E_HP7.Q=((G207$Outflow*-1)+G76$Inflow+S71$Inflow)
  # }else{
  #   L61E_HP7.Q=C41H78$Inflow-((G207$Outflow*-1)+G76$Inflow+S71$Inflow)  
  # }
  L61E_HP7.Q=C41H78$Inflow-((G207$Outflow*-1)+G76$Inflow+S71$Inflow)  
  L61E_HP7.Q=ifelse(L61E_HP7.Q<0,0,L61E_HP7.Q)
  final.dat=data.frame(STRUCT="L61E_HP7",ALIAS="L61E_HP7",Basin="N",WQSite="L61E_HP7",WY=WY(mon.seq$Date.EST[i]),CY=mon.seq$CY[i],month=mon.seq$month[i],Inflow=ifelse(length(L61E_HP7.Q)==0,0,L61E_HP7.Q),Outflow=0)
  L61E_HP7=rbind(L61E_HP7,final.dat)
  print(i)
}
plot(L61E_HP7$Inflow)
flow.mon.sum2=rbind(subset(flow.mon.sum,ALIAS!="C41H78_I"),L61E_HP7)
range(ddply(flow.mon.sum2,"WY",summarise,TFlow.acft=sum(cfs.to.acftd(Inflow),na.rm=T))$TFlow.acft,na.rm=T)

cast(subset(flow.mon.sum2,WY==2018),WY+CY+month~ALIAS,value="Inflow",fun.aggregate=function(x) sum(cfs.to.acftd(x),na.rm=T))

lakeO.WY=ddply(flow.mon.sum2,"WY",summarise,TFlow.in.acft=sum(cfs.to.acftd(Inflow),na.rm=T),TFlow.out.acft=sum(cfs.to.acftd(Outflow),na.rm=T))

lakeO.WY.inflow=cast(flow.mon.sum2,WY~Basin,value="Inflow",fun.aggregate = function(x)sum(cfs.to.acftd(x),na.rm=T))
lakeO.WY.inflow$TFlow=rowSums(lakeO.WY.inflow[,c("E","N","S","W")])
colnames(lakeO.WY.inflow)=c("WY","East","North","South","West","TFlow")

lakeO.WY.outflow=cast(flow.mon.sum2,WY~Basin,value="Outflow",fun.aggregate = function(x)sum(cfs.to.acftd(x),na.rm=T))
lakeO.WY.outflow$TFlow=rowSums(lakeO.WY.outflow[,c("E","N","S","W")])
colnames(lakeO.WY.outflow)=c("WY","East","North","South","West","TFlow")

# Lake Stage
# da.dbks=data.frame(SITE=c("L001","L005","L006","LZ40","S133TW","S352HW","S4TW"),
#                    DBKEY=c("16022","12509","12519","16265","15826","FF579","15732"),type="Daily")

lake.stg=DBHYDRO_daily(dates[1],dates[2],"N3466")
plot(Data.Value~Date,lake.stg)

lake.stg$Date.EST=date.fun(lake.stg$Date)
lake.stg$WY=WY(lake.stg$Date.EST)

stg.dat.WY.sum=ddply(lake.stg,"WY",summarise,mean.val=mean(Data.Value,na.rm=T),min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))

ylim.val=c(9.5,17.5);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(2010,2020);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)

plot(mean.val~WY,stg.dat.WY.sum,ylim=ylim.val)
with(stg.dat.WY.sum,arrows(WY,min.val,WY,max.val,col="grey",angle=90,code=3,length=0.02,lwd=2))

#Lake O rainfall  
lakeO.rf.dbkey=openxlsx::read.xlsx(paste0(data.path,"LakeO_rainfall_dbkey.xlsx"))
paste(lakeO.rf.dbkey$DBKEY,collapse="/")

lakeO.rf=read.csv(paste0(data.path,"20200603_LakeORF.csv"),skip=length(lakeO.rf.dbkey$DBKEY)+2)
lakeO.rf$Date=with(lakeO.rf,as.POSIXct(as.character(Daily.Date),format="%d-%b-%Y",tz="America/New_York"))
lakeO.rf$Date.EST=date.fun(lakeO.rf$Date)
lakeO.rf$WY=WY(lakeO.rf$Date.EST)
lakeO.rf=subset(lakeO.rf,is.na(Date)==F)
lakeO.rf$Data.Value=with(lakeO.rf,ifelse(Data.Value<0,NA,Data.Value))
lakeO.rf=subset(lakeO.rf,Station!="S308_R");#removing questionable data

plot(Data.Value~Date.EST,lakeO.rf)
ddply(lakeO.rf,c("Station"),summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))
ddply(subset(lakeO.rf,Station=="S308_R"),c("Station","Qualifer"),summarise,min.val=min(Data.Value,na.rm=T),max.val=max(Data.Value,na.rm=T))

lakeO.rf.da=ddply(lakeO.rf,c("Date.EST","WY"),summarise,mean.val=mean(Data.Value,na.rm=T),N.val=N(Data.Value))
#plot(N.val~Date.EST,lakeO.rf.da);#sanity check
lakeO.rf.WY=ddply(lakeO.rf.da,"WY",summarise,TRF.in=sum(mean.val,na.rm=T))
lakeO.rf.WY

barplot(lakeO.rf.WY$TRF.in,ylim=c(60,0),space=0,col="skyblue1")
plot(TRF.in~WY,lakeO.rf.WY,ylim=c(60,0))
with(lakeO.rf.WY,segments(WY,rep(0,N(WY)),WY,TRF.in,col="skyblue1",lwd=4))
with(lakeO.rf.WY,points(WY,TRF.in,pch=21,bg="skyblue1",lwd=0.01))

# EAA South
EAA.out.dbkey=data.frame(DBKEY=c("91688","15041","15037","AN602","91012","90973","90934"),
                         SITE=c("S8","S150","S7","S7","G338","G310","G251"))
                         #Inflow_direction=c(1,1,1,1,1,1,1))

# DBK.val=paste("",EAA.out.dbkey$DBKEY,"",collapse="/",sep="")
# SDATE=format(dates[1],"%Y%m%d")
# EDATE=format(dates[2],"%Y%m%d")
# link=paste("http://my.sfwmd.gov/dbhydroplsql/web_io.report_process?v_period=uspec&v_start_date=",SDATE,"&v_end_date=",EDATE,"&v_report_type=format6&v_target_code=file_csv&v_run_mode=onLine&v_js_flag=Y&v_dbkey=",DBK.val,sep="")
# download.file(link,paste0(data.path,"20200603_EAAQ.csv"))
# 
# eaa.q=data.frame()
# for(i in 1:nrow(EAA.out.dbkey)){
#   tmp=DBHYDRO_daily(dates[1],dates[2],EAA.out.dbkey$DBKEY[i])
#   tmp$DBKEY=as.character(EAA.out.dbkey$DBKEY[i])
#   eaa.q=rbind(tmp,eaa.q)
#   print(i)
# }

eaa.q=read.csv(paste0(data.path,"20200603_EAAQ.csv"),skip=length(EAA.out.dbkey$DBKEY)+2)
eaa.q$Date=with(eaa.q,as.POSIXct(as.character(Daily.Date),format="%d-%b-%Y",tz="America/New_York"))
eaa.q$Date.EST=date.fun(eaa.q$Date)
eaa.q$WY=WY(eaa.q$Date.EST)
eaa.q=subset(eaa.q,is.na(Date)==F)
eaa.q=merge(eaa.q,EAA.out.dbkey,"DBKEY")
eaa.q=subset(eaa.q,!(SITE%in%c("G310","G251")))
eaa.q$Data.Value=with(eaa.q,ifelse(Data.Value<0,0,Data.Value))

range(eaa.q$Data.Value,na.rm=T)
ddply(eaa.q,"SITE",summarise,min.val=min(Data.Value,na.rm=T))

eaa.q.WY=ddply(eaa.q,c("WY"),summarise,TFlow.acft=sum(cfs.to.acftd(Data.Value),na.rm=T))
plot(TFlow.acft~WY,eaa.q.WY)


## ENP SRS
SRS.dbkey=data.frame(DBKEY=c("01313","00610","00621","01310","15042","FB752"),
                     SITE=c(paste0("S12",LETTERS[1:4]),"S333","S334"))
srs.q=data.frame()
for(i in 1:nrow(SRS.dbkey)){
  tmp=DBHYDRO_daily(dates[1],dates[2],SRS.dbkey$DBKEY[i])
  tmp$DBKEY=as.character(SRS.dbkey$DBKEY[i])
  srs.q=rbind(tmp,srs.q)
  print(i)
}
srs.q=merge(srs.q,SRS.dbkey,"DBKEY")
srs.q.da.xtab=cast(srs.q,Date~SITE,value="Data.Value",fun.aggregate=function(x) mean(cfs.to.acftd(x),na.rm=T))
srs.q.da.xtab[is.na(srs.q.da.xtab)]=0
srs.q.da.xtab$S333R=with(srs.q.da.xtab,ifelse(S333-S334<0,0,S333-S334))
srs.q.da.xtab$TFlow=rowSums(srs.q.da.xtab[,c("S12A","S12B","S12C","S12D","S333R")],na.rm=T)
srs.q.da.xtab$Date.EST=date.fun(srs.q.da.xtab$Date)
srs.q.da.xtab$WY=WY(srs.q.da.xtab$Date.EST)

srs.q.WY=ddply(srs.q.da.xtab,"WY",summarise,TFlow.acft=sum(TFlow,na.rm=T))
plot(TFlow.acft~WY,srs.q.WY)



# Plot --------------------------------------------------------------------

ylim.val=c(-400,400)*1e4;by.y=200*1e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(2010,2020);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)

#tiff(filename=paste0(plot.path,"LakeO_Discharge.tiff"),width=6.5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,2,0.25,0.5),oma=c(0.25,2.25,0.5,0.75));
layout(matrix(1:2,2,1),heights=c(1,0.25))

cols=c("dodgerblue1","forestgreen","indianred1","grey")#wesanderson::wes_palette("Zissou1",4,"continuous")
plot(TFlow~WY,lakeO.WY.outflow,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA,xaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(lakeO.WY.inflow,shaded.range(WY,rep(0,N(WY)),North,cols[1],lty=1))
with(lakeO.WY.inflow,shaded.range(WY,North,North+East,cols[2],lty=1))
with(lakeO.WY.inflow,shaded.range(WY,North+East,North+East+West,cols[3],lty=1))
with(lakeO.WY.inflow,shaded.range(WY,North+East+West,North+East+West+South,cols[4],lty=1))
with(lakeO.WY.outflow,shaded.range(WY,rep(0,N(WY)),North*-1,cols[1],lty=1))
with(lakeO.WY.outflow,shaded.range(WY,North*-1,(North+East)*-1,cols[2],lty=1))
with(lakeO.WY.outflow,shaded.range(WY,(North+East)*-1,(North+East+West)*-1,cols[3],lty=1))
with(lakeO.WY.outflow,shaded.range(WY,(North+East+West)*-1,(North+East+West+South)*-1,cols[4],lty=1))
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e4)
box(lwd=1)
text(2020.15,-200e4,"Outflow",srt=90,xpd=NA)
text(2020.15,200e4,"Inflow",srt=90,xpd=NA)
mtext(side=2,line=2.5,"Total Discharge (x10\u00B3 Ac-Ft Yr\u207B\u00B9)")
mtext(side=1,line=2,"Florida Water Year")
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
leg.text=c("North","East","West","South")
legend(0.5,0,legend=leg.text,pch=22,pt.bg=adjustcolor(cols,0.25),
       lty=NA,lwd=0.5,col=cols,
       pt.cex=1.5,cex=0.7,ncol=4,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

#tiff(filename=paste0(plot.path,"LakeO_Discharge_bars.tiff"),width=6.5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,2,0.25,0.5),oma=c(0.25,2.25,0.5,0.75));
layout(matrix(1:2,2,1),heights=c(1,0.25))

x=barplot(t(lakeO.WY.inflow[,c(3,2,5,4)]),col=adjustcolor(cols,0.5),space=0,ylim=ylim.val,axes=F,xlab=NA,axisname=F)
abline(h=0)
barplot(t(lakeO.WY.outflow[,c(3,2,5,4)]*-1),col=adjustcolor(cols,0.5),space=0,ylim=ylim.val,add=T,axes=F,xlab=NA,axisname=F)
axis_fun(1,x[seq(1,length(x),2)],x,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e4);box(lwd=1)
mtext(side=2,line=2.5,"Total Discharge (x10\u00B3 Ac-Ft Yr\u207B\u00B9)")
mtext(side=1,line=2,"Florida Water Year")
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
leg.text=c("North","East","West","South")
legend(0.5,0,legend=leg.text,pch=22,pt.bg=adjustcolor(cols,0.5),
       lty=NA,lwd=0.1,col="black",
       pt.cex=1.5,cex=1,ncol=4,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()



#tiff(filename=paste0(plot.path,"LakeO_all.tiff"),width=6.5,height=6.25,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.25,2,0.25,0.25),oma=c(2,4,0.5,0.75),lwd=0.1);
layout(matrix(1:10,5,2,byrow=T),widths=c(1,0.5))

xlim.val=c(2010,2020);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(60,0);by.y=20;ymaj=seq(ylim.val[2],ylim.val[1],by.y);ymin=seq(ylim.val[2],ylim.val[1],by.y/2)
x=barplot(lakeO.rf.WY$TRF.in,ylim=ylim.val,col=adjustcolor("dodgerblue1",0.5),space=0,axes=F,xlab=NA,axisname=F)
axis_fun(1,x[seq(1,length(x),2)],x,NA,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"Annual\nRainfall (in)")
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend="Total average rainfall\nacross Lake Okeechobee",
       pch=22,pt.bg=adjustcolor("dodgerblue1",0.5),
       lty=NA,lwd=0.1,col="black",
       pt.cex=3,cex=1,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(9.5,17.5);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(stg.dat.WY.sum$mean.val,ylim=ylim.val,col=NA,border=NA,space=0,axes=F,xlab=NA,axisname=F)
with(stg.dat.WY.sum,arrows(x,min.val,x,max.val,col="grey",angle=90,code=3,length=0.02,lwd=2))
with(stg.dat.WY.sum,points(x,mean.val,pch=21,bg="grey50"))
axis_fun(1,x[seq(1,length(x),2)],x,NA,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"Mean Stage\n(Ft, NGVD29)")
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend="Annual Average Lake\nStage Elevation\n(\u00B1 min and max)",
       pch=21,pt.bg="grey50",
       lty=NA,lwd=0.1,col="black",
       pt.cex=3,cex=1,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(-400,400)*1e4;by.y=200*1e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(t(lakeO.WY.inflow[,c(3,2,5,4)]),col=adjustcolor(cols,0.5),space=0,ylim=ylim.val,axes=F,xlab=NA,axisname=F)
abline(h=0)
barplot(t(lakeO.WY.outflow[,c(3,2,5,4)]*-1),col=adjustcolor(cols,0.5),space=0,ylim=ylim.val,add=T,axes=F,xlab=NA,axisname=F)
axis_fun(1,x[seq(1,length(x),2)],x,NA,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e4);box(lwd=1)
#mtext(side=2,line=2.5,"Total Discharge (x10\u00B3 Ac-Ft Yr\u207B\u00B9)")
text(x[length(x)]+x[length(x)]*0.1,-200e4,"Outflow",srt=90,xpd=NA)
text(x[length(x)]+x[length(x)]*0.1,200e4,"Inflow",srt=90,xpd=NA)
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
leg.text=c("North","East","West","South")
legend(0.5,0.5,legend=leg.text,pch=22,pt.bg=adjustcolor(cols,0.5),
       lty=NA,lwd=0.1,col="black",
       pt.cex=1.5,cex=1,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,title.adj = 0,title="Lake Okeechobee Regions")

ylim.val=c(0,100)*1e4;by.y=50*1e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(eaa.q.WY$TFlow.acft,col=adjustcolor("grey",0.5),space=0,ylim=ylim.val,axes=F,xlab=NA,axisname=F)
axis_fun(1,x[seq(1,length(x),2)],x,NA,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e4);box(lwd=1)
mtext(side=2,line=2.75,"Total Discharge (x10\u00B3 Ac-Ft Yr\u207B\u00B9)")
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend="Annual Discharge from\ncenteral EAA\n(S-8, S-150, S-7 & G-338)",
       pch=22,pt.bg=adjustcolor("grey",0.5),
       lty=NA,lwd=0.1,col="black",
       pt.cex=3,cex=1,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

ylim.val=c(0,200)*1e4;by.y=100*1e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(srs.q.WY$TFlow.acft,col=adjustcolor("grey",0.5),space=0,ylim=ylim.val,axes=F,xlab=NA,axisname=F)
axis_fun(1,x[seq(1,length(x),2)],x,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj/1e4);box(lwd=1)
mtext(side=1,line=2,"Florida Water Year")
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend="Annual Discharge to ENP\nShark River Slough\n(S12s + S333 - S334)",
       pch=22,pt.bg=adjustcolor("grey",0.5),
       lty=NA,lwd=0.1,col="black",
       pt.cex=3,cex=1,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()