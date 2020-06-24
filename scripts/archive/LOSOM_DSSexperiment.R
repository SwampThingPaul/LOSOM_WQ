# devtools::install_github("eheisman/dssrip",ref="tidyup", INSTALL_opts = "--no-multiarch")
# 
# library(dssrip);

wd="C:/Julian_LaCie/_Github/LOSOM_WQ"

dss_out=opendss(paste0(wd,"/Data/RSMBN/pro_19652016_LORS08/RSMBN_output.dss"))

paths=paste0("/RSMBN/LOK/STAGE/01JAN1965 - 01JAN2016/1DAY/SIMULATED/")
LOK.dat=data.frame(getFullTSC(dss_out,paths))
LOK.dat$Date=date.fun(rownames(LOK.dat))

plot(STAGE~Date,LOK.dat)                      


paths=paste0("/RSMBN/S77/FLOW/01JAN1965 - 01JAN2016/1DAY/SIMULATED/")
S77.dat=data.frame(getFullTSC(dss_out,paths))
S77.dat$Date=date.fun(rownames(S77.dat))

plot(FLOW~Date,S77.dat,type="l")
