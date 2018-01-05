##Script to calculate and write the standard set of derived variables
##for the 800m data

##source('/home/ssobie/assessments/new.netcdf.calendar.R')
##source('/home/ssobie/assessments/cedar.derived.variable.support.R')

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/home/ssobie/code/repos/derived_gridded_indices/cedar.derived.variable.support.r')

library(ncdf4)
library(PCICt)
library(climdex.pcic)

library(doParallel)
registerDoParallel(cores=2)


degree.days.for.model <- function(degree.names,dd.ncs,lat.bnds,n.lat,
                                  tas.subset,tas.list,tasmax.dates) {

    ##Degree Day -could be replaced by the function at the beginning
    yearly.fac <- as.factor(format(tasmax.dates,'%Y'))
    for (k in seq_along(degree.names)) {
      ##print(degree.names[k])
      flag <- is.na(tas.subset[,1])
      degree.name <- degree.names[k]
      dd.nc <- dd.ncs[[k]]
      fx <- dd.fxns[[degree.names[k]]]
      degree.values <- foreach(
                         tas=tas.list,
                         .export=c('yearly.fac','fx')
                         ) %dopar% {
                              degree.values <- fx(tas,yearly.fac)
                         }
      ncol <- length(degree.values[[1]])
      degree.matrix <- matrix(unlist(degree.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
      degree.matrix[flag,] <- NA

      ncvar_put(dd.ncs[[k]],varid=degree.names[k],vals=degree.matrix,
                start=c(i,lat.bnds[1],1),count=c(1,lat.bnds[2],-1))
    }##Degree day loop
}

##----------------------------------------------------------------------------------------------
annual.averages.for.model <- function(ann.names,ann.ncs,lat.bnds,n.lat,
                                      tasmax.subset,tasmax.list,tasmax.dates,
                                      tasmin.subset,tasmin.list,tasmin.dates,
                                      pr.subset,pr.list,pr.dates) {
      
   ##Variables
   for (k in seq_along(ann.names)) {          
       ann.name <- ann.names[k]
       ##print(ann.names[k])
       ann.subset <- switch(ann.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       ann.list <- switch(ann.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)
       ann.dates <- switch(ann.name,tasmax=tasmax.dates,tasmin=tasmin.dates,pr=pr.dates)
       yearly.fac <- as.factor(format(ann.dates,'%Y'))
       flag <- is.na(ann.subset[,1])
       ann.nc <- ann.ncs[[k]]
       ann.fx <- ann.fxns[[ann.name]]
       ann.avg.values <- foreach(
                         data=ann.list,
                         .export=c('yearly.fac','ann.fx')
                         ) %dopar% {
                              ann.avg.values <- ann.fx(data,yearly.fac)
                         }
      ncol <- length(ann.avg.values[[1]])
      ann.avg.matrix <- matrix(unlist(ann.avg.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
      ann.avg.matrix[flag,] <- NA

      ncvar_put(ann.nc,varid=ann.name,vals=ann.avg.matrix,
                start=c(i,lat.bnds[1],1),count=c(1,lat.bnds[2],-1))
    }##Annual average vars loop
}

##----------------------------------------------------------------------------------------------
seasonal.averages.for.model <- function(seas.names,seas.ncs,lat.bnds,n.lat,
                                      tasmax.subset,tasmax.list,tasmax.dates,
                                      tasmin.subset,tasmin.list,tasmin.dates,
                                      pr.subset,pr.list,pr.dates) {
   seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
      
   ##Variables
   for (k in seq_along(seas.names)) {          
       seas.name <- seas.names[k]
       ##print(seas.names[k])
       seas.subset <- switch(seas.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       seas.list <- switch(seas.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)
       seas.dates <- switch(seas.name,tasmax=tasmax.dates,tasmin=tasmin.dates,pr=pr.dates)
       ##Roll over the December months to compute proper seasons
       years <- as.numeric(format(seas.dates,'%Y'))
       uni.yrs <- unique(years)
       months <- as.numeric(format(seas.dates,'%m'))
       dec.ix <- grep(12,months)
       years[dec.ix] <- years[dec.ix] + 1
       dec.fix <- years %in% uni.yrs
       yearly.fac <- as.factor(years[dec.fix])
       monthly.fac <- as.factor(format(seas.dates[dec.fix],'%Y-%m'))      
       seasonal.fac <- factor(seasons[monthly.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
       
       flag <- is.na(seas.subset[,1])
       seas.nc <- seas.ncs[[k]]
       seas.fx <- seas.fxns[[seas.name]]


       avg.fac <- list(yearly.fac,seasonal.fac)       
       seas.corr.list <- lapply(seas.list,function(x,y){x[y]},dec.fix)
       seas.avg.values <- foreach(
                         data=seas.corr.list,
                         .export=c('avg.fac','seas.fx')
                         ) %dopar% {
                              seas.avg.values <- seas.fx(data,avg.fac)
                         }
      ncol <- length(seas.avg.values[[1]])
      seas.avg.matrix <- matrix(unlist(seas.avg.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
      seas.avg.matrix[flag,] <- NA
      ncvar_put(seas.nc,varid=seas.name,vals=seas.avg.matrix,
                start=c(i,lat.bnds[1],1),count=c(1,lat.bnds[2],-1))
    }##Seasonal average vars loop
}

##----------------------------------------------------------------------------------------------
monthly.averages.for.model <- function(mon.names,mon.ncs,lat.bnds,n.lat,
                                      tasmax.subset,tasmax.list,tasmax.dates,
                                      tasmin.subset,tasmin.list,tasmin.dates,
                                      pr.subset,pr.list,pr.dates) {
      
   ##Variables
   for (k in seq_along(mon.names)) {          
       mon.name <- mon.names[k]
       ##print(mon.names[k])
       mon.subset <- switch(mon.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       mon.list <- switch(mon.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)
       mon.dates <- switch(mon.name,tasmax=tasmax.dates,tasmin=tasmin.dates,pr=pr.dates)
       monthly.fac <- as.factor(format(mon.dates,'%Y-%m'))
       flag <- is.na(mon.subset[,1])
       mon.nc <- mon.ncs[[k]]
       mon.fx <- mon.fxns[[mon.name]]
       mon.avg.values <- foreach(
                         data=mon.list,
                         .export=c('monthly.fac','mon.fx')
                         ) %dopar% {
                              mon.avg.values <- mon.fx(data,monthly.fac)
                         }
      ncol <- length(mon.avg.values[[1]])
      mon.avg.matrix <- matrix(unlist(mon.avg.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
      mon.avg.matrix[flag,] <- NA

      ncvar_put(mon.nc,varid=mon.name,vals=mon.avg.matrix,
                start=c(i,lat.bnds[1],1),count=c(1,lat.bnds[2],-1))
    }##Monthly average vars loop
}

##----------------------------------------------------------------------------------------------
annual.extremes.for.model <- function(ext.names,ext.ncs,lat.bnds,n.lat,
                                      tasmax.subset,tasmax.list,tasmax.dates,
                                      tasmin.subset,tasmin.list,tasmin.dates,
                                      pr.subset,pr.list,pr.dates) {
      
   ##Variables
   for (k in seq_along(ext.names)) {          
       ext.name <- ext.names[k]
       ##print(ext.names[k])
       ext.subset <- switch(ext.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       ext.list <- switch(ext.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)
       ext.dates <- switch(ext.name,tasmax=tasmax.dates,tasmin=tasmin.dates,pr=pr.dates)
       yearly.fac <- as.factor(format(ext.dates,'%Y'))
       flag <- is.na(ext.subset[,1])
       ext.nc <- ext.ncs[[k]]
       ext.fx <- ext.fxns[[ext.name]]
       ext.avg.values <- foreach(
                         data=ext.list,
                         .export=c('yearly.fac','ext.fx')
                         ) %dopar% {
                              ext.avg.values <- ext.fx(data,yearly.fac)
                         }
      ncol <- length(ext.avg.values[[1]])
      ext.avg.matrix <- matrix(unlist(ext.avg.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
      ext.avg.matrix[flag,] <- NA

      ncvar_put(ext.nc,varid=ext.name,vals=ext.avg.matrix,
                start=c(i,lat.bnds[1],1),count=c(1,lat.bnds[2],-1))
    }##Annual extremes vars loop
}

##----------------------------------------------------------------------------------------------
climdex.extremes.for.model <- function(clim.names,clim.ncs,lat.bnds,n.lat,
                                      tasmax.subset,tasmax.list,tasmax.dates,
                                      tasmin.subset,tasmin.list,tasmin.dates,
                                      pr.subset,pr.list,pr.dates,base) {

   flag <- is.na(tasmax.subset[,1])      
   climdex.objects <- foreach(
                             tasmax=tasmax.list,
                             tasmin=tasmin.list,
                             pr=pr.list,
                             .export=c('climdexInput.raw','tasmax.dates','tasmin.dates','pr.dates','base')
                             ) %dopar% {
                                objects <- climdexInput.raw(tmax=tasmax,tmin=tasmin,prec=pr,
                                           tmax.dates=tasmax.dates,tmin.dates=tasmin.dates,prec.dates=pr.dates,base=base)
                             }

   for (k in seq_along(clim.names)) {
        clim.name <- clim.names[k]        
        ##print(clim.names[k])
        clim.nc <- clim.ncs[[k]]
        clim.fx <- clim.fxns[[clim.name]]
        climdex.info <- get.climdex.info(clim.name)
        climdex.values <- foreach(
                                 obj=climdex.objects,
                                 .export='clim.fx'
                                 ) %dopar% {
                                 climdex.values <- clim.fx(obj)
                                 }
        ncol <- length(climdex.values[[1]])
        climdex.matrix <- matrix(unlist(climdex.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
        climdex.matrix[flag,] <- NA
        ncvar_put(clim.nc,varid=climdex.info[1],vals=climdex.matrix,
                  start=c(i,lat.bnds[1],1),count=c(1,lat.bnds[2],-1))
    }##Climdex names vars loop

}





##--------------------------------------------------------------


##****************************************************************

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}
tmp.dir <- '/local_temp/ssobie/prism/' ##tmpdir


gcm <- 'CanESM2'
scenario <- 'rcp85'
run <- 'r1i1p1'
interval <- '1951-1954'

##Latitude Bands
lat.st <- format(seq(48.0,59.9,0.1),nsmall=1)
lat.en <- format(seq(48.1,60.0,0.1),nsmall=1)

len <- length(lat.st)

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/prism/'
data.dir <- paste0(base.dir,gcm,'/lat_split/')
derived.dir <- paste0(base.dir,gcm,'/',scenario,'/')

##Move data to local storage for better I/O
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

##---------------------------------------------------------------------------
##Degree Day Files for writing
degree.names <- c('cdd','fdd','gdd','hdd')
dd.dir <- paste0(base.dir,gcm,'/',scenario,'/degree_days/')
dd.ncs <- vector(mode='list',length=length(degree.names))
for (d in seq_along(degree.names)) {
  dd.file <- paste0(dd.dir,degree.names[d],'_annual_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  dd.ncs[[d]] <- nc_open(dd.file,write=TRUE)
}
common.lat <- ncvar_get(dd.ncs[[1]],'lat')

##---------------------------------------------------------------------------
##Annual Average Files for writing
ann.names <- c('pr','tasmax','tasmin')
avg.type <- c('total','average','average')
ann.dir <- paste0(base.dir,gcm,'/',scenario,'/annual/')
ann.ncs <- vector(mode='list',length=length(ann.names))
for (a in seq_along(ann.names)) {
  ann.file <- paste0(ann.dir,ann.names[a],'_annual_',avg.type[a],'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  ann.ncs[[a]] <- nc_open(ann.file,write=TRUE)
 
}
common.lat <- ncvar_get(ann.ncs[[1]],'lat')

##---------------------------------------------------------------------------
##Annual Block Maxima Files for writing
ext.names <- c('pr','tasmax','tasmin')
ext.type <- c('maximum','maximum','minimum')
ext.dir <- paste0(base.dir,gcm,'/',scenario,'/annual_extremes/')
ext.ncs <- vector(mode='list',length=length(ext.names))
for (a in seq_along(ext.names)) {
  ext.file <- paste0(ext.dir,ext.names[a],'_annual_',ext.type[a],'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  ext.ncs[[a]] <- nc_open(ext.file,write=TRUE)
 
}
common.lat <- ncvar_get(ext.ncs[[1]],'lat')


##---------------------------------------------------------------------------
##Seasonal Average Files for writing
seas.names <- c('pr','tasmax','tasmin')
seas.dir <- paste0(base.dir,gcm,'/',scenario,'/seasonal/')
seas.ncs <- vector(mode='list',length=length(seas.names))
for (a in seq_along(seas.names)) {
  seas.file <- paste0(seas.dir,seas.names[a],'_seasonal_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  seas.ncs[[a]] <- nc_open(seas.file,write=TRUE)
 
}
common.lat <- ncvar_get(seas.ncs[[1]],'lat')

##---------------------------------------------------------------------------
##Monthly Average Files for writing
mon.names <- c('pr','tasmax','tasmin')
mon.dir <- paste0(base.dir,gcm,'/',scenario,'/monthly/')
mon.ncs <- vector(mode='list',length=length(mon.names))
for (a in seq_along(mon.names)) {
  mon.file <- paste0(mon.dir,mon.names[a],'_monthly_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  mon.ncs[[a]] <- nc_open(mon.file,write=TRUE)
 
}
common.lat <- ncvar_get(mon.ncs[[1]],'lat')

##---------------------------------------------------------------------------
##Climdex Files for writing
clim.names <- c('climdex.fd','climdex.su','climdex.id',
                'climdex.tr','climdex.gsl','climdex.txx',
                'climdex.tnx','climdex.txn','climdex.tnn',
                'climdex.tn10p','climdex.tx10p','climdex.tn90p',
                'climdex.tx90p','climdex.wsdi','climdex.csdi',
                'climdex.dtr','climdex.rx1day','climdex.rx2day',
                'climdex.rx5day','climdex.sdii','climdex.r10mm',
                'climdex.r20mm','climdex.cdd','climdex.cwd',
                'climdex.r95ptot','climdex.r99ptot','climdex.prcptot',
                'climdex.r95days','climdex.r99days','climdex.r95store',
                'climdex.r99store')

climdex.dir <- paste0(base.dir,gcm,'/',scenario,'/climdex/')
clim.ncs <- vector(mode='list',length=length(clim.names))
for (d in seq_along(clim.names)) {
  clim.name <- clim.names[d]
  climdex.info <- get.climdex.info(clim.name)
  climdex.var <- climdex.info[1]
  climdex.calendar <- climdex.info[2]
  clim.file <- paste0(climdex.dir,climdex.var,'_',tolower(climdex.calendar),'_BCCAQ2-PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  clim.ncs[[d]] <- nc_open(clim.file,write=TRUE)
}
common.lat <- ncvar_get(clim.ncs[[1]],'lat')

##---------------------------------------------------------------------------
##---------------------------------------------------------------------------

##Iterate over the latitude files
for (i in 1:2) { ##len) {
  lat.interval <- paste0(lat.st[i],'-',lat.en[i])
  tasmax.file <- paste0('tasmax_gcm_prism_BCCAQ_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmax.nc <- nc_open(paste0(data.dir,tasmax.file),write=FALSE)
  tasmax.dates <- netcdf.calendar(tasmax.nc)

  tasmin.file <- paste0('tasmin_gcm_prism_BCCAQ_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmin.nc <- nc_open(paste0(data.dir,tasmin.file),write=FALSE)
  tasmin.dates <- netcdf.calendar(tasmin.nc)

  pr.file <- paste0('pr_gcm_prism_BCCAQ_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  pr.nc <- nc_open(paste0(data.dir,pr.file),write=FALSE)
  pr.dates <- netcdf.calendar(pr.nc)

  lon <- ncvar_get(tasmax.nc,'lon')
  lat <- ncvar_get(tasmax.nc,'lat')
  n.lon <- length(lon)
  n.lat <- length(lat)
  lat.match <- which(common.lat %in% lat)
  lat.bnds <- c(lat.match[1],(tail(lat.match,1)-lat.match[1])+1)
  print('Latitude bands:')
  print(lat.bnds)

  for (i in 1:n.lon) {
    print(paste0('Longitude: ',i,' of ',n.lon))
    tasmax.subset <- ncvar_get(tasmax.nc,'tasmax',start=c(i,1,1),count=c(1,-1,-1))
    tasmin.subset <- ncvar_get(tasmin.nc,'tasmin',start=c(i,1,1),count=c(1,-1,-1))    
    pr.subset <- ncvar_get(pr.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))
    tas.subset <- (tasmax.subset + tasmin.subset)/2

    tasmax.list <- list()
    tasmin.list <- list()
    tas.list <- list()
    pr.list <- list()

    for (j in 1:n.lat) {
      tasmax.list[[j]] <- tasmax.subset[j,]
      tasmin.list[[j]] <- tasmin.subset[j,]
      tas.list[[j]] <- tas.subset[j,]
      pr.list[[j]] <- pr.subset[j,]
    }

    ##----------------------------------------------------------
    ##Degree Day 
      degree.days.for.model(degree.names,dd.ncs,lat.bnds,n.lat,
                            tas.subset,tas.list,tasmax.dates)
    ##----------------------------------------------------------
    ##Annual Averages 
      annual.averages.for.model(ann.names,ann.ncs,lat.bnds,n.lat,
                                tasmax.subset,tasmax.list,tasmax.dates,
                                tasmin.subset,tasmin.list,tasmin.dates,
                                pr.subset,pr.list,pr.dates)
    ##----------------------------------------------------------
    ##Seasonal Averages 
      seasonal.averages.for.model(seas.names,seas.ncs,lat.bnds,n.lat,
                                tasmax.subset,tasmax.list,tasmax.dates,
                                tasmin.subset,tasmin.list,tasmin.dates,
                                pr.subset,pr.list,pr.dates)
    ##----------------------------------------------------------
    ##Monthly Averages 
      monthly.averages.for.model(mon.names,mon.ncs,lat.bnds,n.lat,
                                tasmax.subset,tasmax.list,tasmax.dates,
                                tasmin.subset,tasmin.list,tasmin.dates,
                                pr.subset,pr.list,pr.dates)
    ##----------------------------------------------------------
    ##Annual Extremes
      annual.extremes.for.model(ext.names,ext.ncs,lat.bnds,n.lat,
                                tasmax.subset,tasmax.list,tasmax.dates,
                                tasmin.subset,tasmin.list,tasmin.dates,
                                pr.subset,pr.list,pr.dates)
    ##----------------------------------------------------------
    ##Climdex
      climdex.extremes.for.model(clim.names,clim.ncs,lat.bnds,n.lat,
                                 tasmax.subset,tasmax.list,tasmax.dates,
                                 tasmin.subset,tasmin.list,tasmin.dates,
                                 pr.subset,pr.list,pr.dates,base=c(1951,1954))
    ##----------------------------------------------------------

  }##Longitude Loop
  nc_close(tasmax.nc)
  nc_close(tasmin.nc)
  nc_close(pr.nc)
}##Latitude File Loop

for (d in seq_along(degree.names)) {
  nc_close(dd.ncs[[d]])
}

for (d in seq_along(ann.names)) {
  nc_close(ann.ncs[[d]])
}

for (d in seq_along(seas.names)) {
  nc_close(seas.ncs[[d]])
}

for (d in seq_along(mon.names)) {
  nc_close(mon.ncs[[d]])
}

for (d in seq_along(ext.names)) {
  nc_close(ext.ncs[[d]])
}

for (d in seq_along(clim.names)) {
  nc_close(clim.ncs[[d]])
}
