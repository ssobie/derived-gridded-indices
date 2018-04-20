##Script to calculate and write the standard set of derived variables
##for the 800m data

##Rprof('cedar.derived.out')

ptm <- proc.time()
cedar <- FALSE

if (cedar) {
   source('/home/ssobie/assessments/new.netcdf.calendar.R')
   source('/home/ssobie/assessments/cedar.derived.variable.support.r')
} else {
   source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
   source('/home/ssobie/code/repos/derived_gridded_indices/cedar.derived.variable.support.r')
}

library(ncdf4)
library(PCICt)
library(climdex.pcic)

library(foreach)
library(doParallel)
registerDoParallel(cores=2)

get.seasonal.fac <- function(seas.dates) {
   seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
   years <- as.numeric(format(seas.dates,'%Y'))
   uni.yrs <- unique(years)
   months <- as.numeric(format(seas.dates,'%m'))
   dec.ix <- grep(12,months)
   years[dec.ix] <- years[dec.ix] + 1
   dec.fix <- years %in% uni.yrs
   yearly.fac <- as.factor(years[dec.fix])
   monthly.fac <- as.factor(format(seas.dates[dec.fix],'%m'))      
   seasonal.fac <- factor(seasons[monthly.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
   avg.fac <- list(yearly.fac,seasonal.fac)       

  return(list(fac=avg.fac,fix=dec.fix))
}


degree.days.for.model <- function(degree.names,dd.ncs,lat.ix,n.lon,
                                  tas.subset,tas.list,yearly.fac) {

    ##Degree Day -could be replaced by the function at the beginning
    for (k in seq_along(degree.names)) {
      print(degree.names[k])
      flag <- is.na(tas.subset[,1])
      flen <- sum(!flag)
      print('Number of valid cells:')
      print(flen)
      degree.name <- degree.names[k]
      dd.nc <- dd.ncs[[k]]
      fx <- dd.fxns[[degree.names[k]]]
      sub.list <- tas.list[!flag]

      ftm <- proc.time()      
      sub.values <- foreach(
                         temp=sub.list,
                         .export=c('yearly.fac','fx','dd')
                         ) %dopar% {
                              degree.values <- fx(temp,yearly.fac)
                              degree.values
                         }
       print('Foreach time:')
       print(proc.time()-ftm)

      ncol <- length(sub.values[[1]])
      degree.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(sub.values),nrow=flen,ncol=ncol,byrow=TRUE)
      degree.matrix[!flag,] <- sub.matrix
      rm(sub.values)
      rm(sub.matrix)
      ncvar_put(dd.ncs[[k]],varid=degree.names[k],vals=degree.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
      rm(degree.matrix)       
    }##Degree day loop
}

##----------------------------------------------------------------------------------------------
annual.averages.for.model <- function(ann.names,ann.ncs,lat.ix,n.lon,yearly.fac,
                                      tasmax.subset,tasmax.list,
                                      tasmin.subset,tasmin.list,
                                      pr.subset,pr.list) {
                                            
   ##Variables
   for (k in seq_along(ann.names)) {          
       ann.name <- ann.names[k]
       print(ann.names[k])
       ann.subset <- switch(ann.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       ann.list <- switch(ann.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)
       flag <- is.na(ann.subset[,1])
       flen <- sum(!flag)
       ann.nc <- ann.ncs[[k]]
       ann.fx <- ann.fxns[[ann.name]]
       sub.list <- ann.list[!flag]

       ftm <- proc.time()       	      
       ann.avg.values <- foreach(
                         data=sub.list,
                         .export=c('yearly.fac','ann.fx')
                         ) %dopar% {
                              ann.avg.values <- ann.fx(data,yearly.fac)
                              ann.avg.values
                         }
       ##print('Foreach time:')
       ##print(proc.time()-ftm)	      
      ncol <- length(ann.avg.values[[1]])
      ann.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(ann.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
      rm(ann.avg.values)
      ann.avg.matrix[!flag] <- sub.matrix
      rm(sub.matrix)
      print(dim(ann.avg.matrix))
      ncvar_put(ann.nc,varid=ann.name,vals=ann.avg.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
    }##Annual average vars loop
}

##----------------------------------------------------------------------------------------------
seasonal.averages.for.model <- function(seas.names,seas.ncs,lat.ix,n.lon,seasonal.fac,
                                      tasmax.subset,tasmax.list,
                                      tasmin.subset,tasmin.list,
                                      pr.subset,pr.list) {

   ##Variables
   for (k in seq_along(seas.names)) {          
       seas.name <- seas.names[k]
       print(seas.names[k])
       seas.subset <- switch(seas.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       seas.list <- switch(seas.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)

       ##Roll over the December months to compute proper seasons
       flag <- is.na(seas.subset[,1])
       flen <- sum(!flag)
       seas.nc <- seas.ncs[[k]]
       seas.fx <- seas.fxns[[seas.name]]
       avg.fac <- seasonal.fac$fac
       dec.fix <- seasonal.fac$fix
       seas.corr.list <- lapply(seas.list,function(x,y){x[y]},dec.fix)
       sub.list <- seas.corr.list[!flag]

       ftm <- proc.time()       
       seas.avg.values <- foreach(
                         data=sub.list,
                         .export=c('avg.fac','seas.fx')
                         ) %dopar% {
                              seas.avg.values <- as.vector(t(seas.fx(data,avg.fac)))
                         }
       print('Foreach time:')
       print(proc.time()-ftm)
			 
      ncol <- length(seas.avg.values[[1]])
      seas.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(seas.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
      seas.avg.matrix[!flag,] <- sub.matrix
      rm(seas.avg.values)
      rm(sub.matrix)
      ncvar_put(seas.nc,varid=seas.name,vals=seas.avg.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
    }##Seasonal average vars loop
}

##----------------------------------------------------------------------------------------------
monthly.averages.for.model <- function(mon.names,mon.ncs,lat.ix,n.lon,monthly.fac,
                                      tasmax.subset,tasmax.list,
                                      tasmin.subset,tasmin.list,
                                      pr.subset,pr.list) {
      
   ##Variables
   for (k in seq_along(mon.names)) {          
       mon.name <- mon.names[k]
       print(mon.names[k])
       mon.subset <- switch(mon.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       mon.list <- switch(mon.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)

       flag <- is.na(mon.subset[,1])
       flen <- sum(!flag)
       mon.nc <- mon.ncs[[k]]
       mon.fx <- mon.fxns[[mon.name]]
       sub.list <- mon.list[!flag]

       ftm <- proc.time()       
       mon.avg.values <- foreach(
                         data=sub.list,
                         .export=c('monthly.fac','mon.fx')
                         ) %dopar% {
                              mon.avg.values <- mon.fx(data,monthly.fac)
                         }
       print('Foreach time:')
       print(proc.time()-ftm)
			 
      ncol <- length(mon.avg.values[[1]])
      mon.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)      
      sub.matrix <- matrix(unlist(mon.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
      mon.avg.matrix[!flag,] <- sub.matrix
      rm(mon.avg.values)
      rm(sub.matrix)
      ncvar_put(mon.nc,varid=mon.name,vals=mon.avg.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
    }##Monthly average vars loop
}

##----------------------------------------------------------------------------------------------
annual.extremes.for.model <- function(ext.names,ext.ncs,lat.ix,n.lon,yearly.fac,
                                      tasmax.subset,tasmax.list,
                                      tasmin.subset,tasmin.list,
                                      pr.subset,pr.list) {
      
   ##Variables
   for (k in seq_along(ext.names)) {          
       ext.name <- ext.names[k]
       print(ext.names[k])
       ext.subset <- switch(ext.name,tasmax=tasmax.subset,tasmin=tasmin.subset,pr=pr.subset)
       ext.list <- switch(ext.name,tasmax=tasmax.list,tasmin=tasmin.list,pr=pr.list)
       flag <- is.na(ext.subset[,1])
       flen <- sum(!flag)
       sub.list <- ext.list[!flag]

       ext.nc <- ext.ncs[[k]]
       ext.fx <- ext.fxns[[ext.name]]

       ftm <- proc.time()       
       ext.avg.values <- foreach(
                         data=sub.list,
                         .export=c('yearly.fac','ext.fx')
                         ) %dopar% {
                              ext.avg.values <- ext.fx(data,yearly.fac)
                         }
       print('Foreach time:')
       print(proc.time()-ftm)
			 
      ncol <- length(ext.avg.values[[1]])
      ext.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(ext.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
      ext.avg.matrix[!flag,] <- sub.matrix
      rm(ext.avg.values)
      rm(sub.matrix)
      ncvar_put(ext.nc,varid=ext.name,vals=ext.avg.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
      rm(ext.avg.matrix)		
    }##Annual extremes vars loop
}

##----------------------------------------------------------------------------------------------
climdex.extremes.for.model <- function(clim.names,clim.ncs,lat.ix,n.lon,
                                      tasmax.subset,tasmax.list,tasmax.dates,
                                      tasmin.subset,tasmin.list,tasmin.dates,
                                      pr.subset,pr.list,pr.dates,base) {

   flag <- is.na(tasmax.subset[,1])      
   flen <- sum(!flag)
   tasmax.sub <- tasmax.list[!flag]
   tasmin.sub <- tasmin.list[!flag]
   pr.sub <- pr.list[!flag]

   climdex.objects <- foreach(
                             tasmax=tasmax.sub,
                             tasmin=tasmin.sub,
                             pr=pr.sub,
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
        climdex.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
        sub.matrix <- matrix(unlist(climdex.values),nrow=flen,ncol=ncol,byrow=TRUE)
        rm(climdex.values)
        climdex.matrix[!flag,] <- sub.matrix
        rm(sub.matrix)
        ncvar_put(clim.nc,varid=climdex.info[1],vals=climdex.matrix,
                  start=c(1,lat.ix,1),count=c(-1,1,-1))
        rm(climdex.matrix)                  
    }##Climdex names vars loop
}





##--------------------------------------------------------------
##****************************************************************

if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

tmp.dir <- tmpdir ##'/local_temp/ssobie/prism/' ##tmpdir

##gcm <- 'CNRM-CM5'
##scenario <- 'rcp85'
##run <- 'r1i1p1'
##interval <- '1951-2100'
##type <- 'climdex'


##Latitude Bands
lat.st <- seq(48.1,59.9,0.1) ##format(seq(48.0,59.9,0.1),nsmall=1)
lat.en <- seq(48.2,60.0,0.1) ##format(seq(48.1,60.0,0.1),nsmall=1)

len <- length(lat.st)

if (cedar) {
   base.dir <- '/scratch/ssobie/prism/'
} else {
   base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
}

data.dir <- paste0(base.dir,gcm,'/lat_split/')
template.dir <- paste0(base.dir,gcm,'/template/',scenario,'/',type,'/')
##Move data to local storage for better I/O
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}
dir.create(paste0(tmp.dir,scenario),recursive=T)

rtm <- proc.time()
print('Copying')

file.copy(from=template.dir,to=paste0(tmp.dir,scenario,'/'),overwrite=TRUE,recursive=TRUE) ##
##test.file <- paste0(base.dir,gcm,'/template/',scenario,'/climdex_mon/txxETCCDI_mon_BCCAQ2-PRISM_CNRM-CM5_rcp85_r1i1p1_1951-2100.nc')
##file.copy(from=test.file,to=paste0(tmp.dir,scenario,'/climdex/'),overwrite=TRUE) ##

move.to <- paste("rsync -av ",template.dir," ",tmp.dir,scenario,sep='')
##print(move.to)
##system(move.to)

print('Move template time')
print(proc.time()-rtm)

##---------------------------------------------------------------------------
  ##Degree Day Files for writing
if (type=='degree_days') {
  print('Degree days opening')
  degree.names <- c('cdd','fdd','gdd','hdd')
  dd.dir <- paste0(tmp.dir,scenario,'/degree_days/')
  dd.ncs <- vector(mode='list',length=length(degree.names))
  for (d in seq_along(degree.names)) {
    dd.file <- paste0(dd.dir,degree.names[d],'_annual_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
    dd.ncs[[d]] <- nc_open(dd.file,write=TRUE)
  }
  common.lat <- ncvar_get(dd.ncs[[1]],'lat')
}
##---------------------------------------------------------------------------
if (type=='annual') {
  ##Annual Average Files for writing
  print('Ann Avg opening')
  ann.names <- c('pr','tasmax','tasmin')
  avg.type <- c('total','average','average')
  ann.dir <- paste0(tmp.dir,scenario,'/annual/')
  ann.ncs <- vector(mode='list',length=length(ann.names))
  for (a in seq_along(ann.names)) {
    ann.file <- paste0(ann.dir,ann.names[a],'_annual_',avg.type[a],'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
    ann.ncs[[a]] <- nc_open(ann.file,write=TRUE) 
  }
  common.lat <- ncvar_get(ann.ncs[[1]],'lat')
}
##---------------------------------------------------------------------------
if (type=='annual_extremes') {
##Annual Block Maxima Files for writing
  print('Ann extremes opening')
  ext.names <- c('pr','tasmax','tasmin')
  ext.type <- c('maximum','maximum','minimum')
  ext.dir <- paste0(tmp.dir,scenario,'/annual_extremes/')
  ext.ncs <- vector(mode='list',length=length(ext.names))
  for (a in seq_along(ext.names)) {
    ext.file <- paste0(ext.dir,ext.names[a],'_annual_',ext.type[a],'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
    ext.ncs[[a]] <- nc_open(ext.file,write=TRUE)
  }
  common.lat <- ncvar_get(ext.ncs[[1]],'lat')
}
##---------------------------------------------------------------------------
if (type=='seasonal') {
##Seasonal Average Files for writing
  print('Seasonal averages opening')
  seas.names <- c('pr','tasmax','tasmin')
  seas.dir <- paste0(tmp.dir,scenario,'/seasonal/')
  seas.ncs <- vector(mode='list',length=length(seas.names))
  for (a in seq_along(seas.names)) {
    seas.file <- paste0(seas.dir,seas.names[a],'_seasonal_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
    seas.ncs[[a]] <- nc_open(seas.file,write=TRUE)
  }
  common.lat <- ncvar_get(seas.ncs[[1]],'lat')
}
##---------------------------------------------------------------------------
if (type=='monthly') {
  ##Monthly Average Files for writing
  print('monthly avg opening')
  mon.names <- c('pr','tasmax','tasmin')
  mon.dir <- paste0(tmp.dir,scenario,'/monthly/')
  mon.ncs <- vector(mode='list',length=length(mon.names))
  for (a in seq_along(mon.names)) {
    mon.file <- paste0(mon.dir,mon.names[a],'_monthly_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
    mon.ncs[[a]] <- nc_open(mon.file,write=TRUE)
  }
  common.lat <- ncvar_get(mon.ncs[[1]],'lat')
}
##---------------------------------------------------------------------------
if (type=='climdex') {
##Climdex Files for writing
  print('Climdex opening')
##'climdex.txx','climdex.tnx','climdex.txn','climdex.tnn',
##'climdex.tn10p','climdex.tx10p','climdex.tn90p','climdex.tx90p',
##'climdex.dtr','climdex.rx1day','climdex.rx2day','climdex.rx5day',

##  clim.names <- c('climdex.fd','climdex.su','climdex.id',
##                  'climdex.tr','climdex.gsl','climdex.wsdi','climdex.csdi',
##                  'climdex.sdii','climdex.r10mm',
##                  'climdex.r20mm','climdex.cdd','climdex.cwd',
##                  'climdex.r95ptot','climdex.r99ptot','climdex.prcptot',
##                  'climdex.r95days','climdex.r99days','climdex.r95store',
##                  'climdex.r99store')

clim.names <- 'climdex.txx'

  climdex.dir <- paste0(tmp.dir,'/',scenario,'/climdex/')
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
}
##---------------------------------------------------------------------------
##---------------------------------------------------------------------------

##Iterate over the latitude files
for (i in 1:len) {
  lat.interval <- paste0(sprintf('%.1f',lat.st[i]),'-',sprintf('%.1f',lat.en[i]))

  tasmax.input <- paste0('tasmax_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  print(paste0('Copy ',tasmax.input))
  file.copy(paste0(data.dir,"/",tasmax.input),tmp.dir,overwrite=TRUE)

  tasmin.input <- paste0('tasmin_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  print(paste0('Copy ',tasmin.input))
  file.copy(paste0(data.dir,"/",tasmin.input),tmp.dir,overwrite=TRUE)

  pr.input <- paste0('pr_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  print(paste0('Copy ',pr.input))
  file.copy(paste0(data.dir,"/",pr.input),tmp.dir,overwrite=TRUE)

  print('TX opening')
  tasmax.file <- paste0('tasmax_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmax.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
  tasmax.dates <- netcdf.calendar(tasmax.nc)
  yearly.fac <- as.factor(format(tasmax.dates,'%Y'))
  seasonal.fac <- get.seasonal.fac(tasmax.dates)
  monthly.fac <- as.factor(format(tasmax.dates,'%Y-%m'))  

  print('TN Opening')
  tasmin.file <- paste0('tasmin_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmin.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)
  tasmin.dates <- netcdf.calendar(tasmin.nc)

  print('PR Opening')
  pr.file <- paste0('pr_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  pr.nc <- nc_open(paste0(tmp.dir,pr.file),write=FALSE)
  pr.dates <- netcdf.calendar(pr.nc)
  
  lon <- ncvar_get(tasmax.nc,'lon')
  lat <- ncvar_get(tasmax.nc,'lat')
  n.lon <- length(lon)
  n.lat <- length(lat)
  lat.match <- which(common.lat %in% lat)
  lat.bnds <- c(lat.match[1],(tail(lat.match,1)-lat.match[1])+1)

  print('Latitude bands:')
  print(lat.bnds)

  for (i in 1:n.lat) { ##n.lon) {
    lat.ix <- lat.match[i]
    ltm <- proc.time()

    print(paste0('Latitude: ',i,' of ',n.lat))
    tasmax.subset <- ncvar_get(tasmax.nc,'tasmax',start=c(1,i,1),count=c(-1,1,-1))
    tasmin.subset <- ncvar_get(tasmin.nc,'tasmin',start=c(1,i,1),count=c(-1,1,-1))    
    pr.subset <- ncvar_get(pr.nc,'pr',start=c(1,i,1),count=c(-1,1,-1))
    tas.subset <- (tasmax.subset + tasmin.subset)/2

    tasmax.list <- vector(mode='list',length=n.lon)
    tasmin.list <- vector(mode='list',length=n.lon)
    tas.list <- vector(mode='list',length=n.lon)
    pr.list <- vector(mode='list',length=n.lon)

    rtm <- proc.time()     
    tasmax.list <- lapply(seq_len(nrow(tasmax.subset)), function(i) tasmax.subset[i,])
    tasmin.list <- lapply(seq_len(nrow(tasmin.subset)), function(i) tasmin.subset[i,])     
    pr.list <- lapply(seq_len(nrow(pr.subset)), function(i) pr.subset[i,])
    tas.list <- lapply(seq_len(nrow(tas.subset)), function(i) tas.subset[i,])
    print('Lapply convert to list')
    print(proc.time()-rtm) 

    ##----------------------------------------------------------
    ##Degree Day 
    if (type=='degree_days') {
    rtm <- proc.time()     
       degree.days.for.model(degree.names,dd.ncs,lat.ix,n.lon,
                              tas.subset,tas.list,yearly.fac)
    print('Degree Days time:')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    ##Annual Averages 
    if (type=='annual') {
    rtm <- proc.time()     
      annual.averages.for.model(ann.names,ann.ncs,lat.ix,n.lon,yearly.fac,
                                tasmax.subset,tasmax.list,
                                tasmin.subset,tasmin.list,
                                pr.subset,pr.list)
    print('Annual Averages Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    ##Seasonal Averages 
    if (type=='seasonal') {
    rtm <- proc.time()     
      seasonal.averages.for.model(seas.names,seas.ncs,lat.ix,n.lon,seasonal.fac,
                                tasmax.subset,tasmax.list,
                                tasmin.subset,tasmin.list,
                                pr.subset,pr.list)
    print('Seasonal Averages Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    ##Monthly Averages 
    if (type=='monthly') {
    rtm <- proc.time()     
      monthly.averages.for.model(mon.names,mon.ncs,lat.ix,n.lon,monthly.fac,
                                tasmax.subset,tasmax.list,
                                tasmin.subset,tasmin.list,
                                pr.subset,pr.list)
    print('Monthly Averages Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    if (type=='annual_extremes') {
    ##Annual Extremes
    rtm <- proc.time()     
      annual.extremes.for.model(ext.names,ext.ncs,lat.ix,n.lon,yearly.fac,
                                tasmax.subset,tasmax.list,
                                tasmin.subset,tasmin.list,
                                pr.subset,pr.list)
    print('Annual Extremes Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    if (type=='climdex') {
    ##Climdex
    rtm <- proc.time()     
      climdex.extremes.for.model(clim.names,clim.ncs,lat.ix,n.lon,
                                 tasmax.subset,tasmax.list,tasmax.dates,
                                 tasmin.subset,tasmin.list,tasmin.dates,
                                 pr.subset,pr.list,pr.dates,base=c(1971,2000))
    print('Climdex Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    rm(tasmax.list)
    rm(tasmin.list)
    rm(pr.list)
    rm(tas.list)
    rm(tas.subset)
    rm(tasmax.subset)
    rm(tasmin.subset)
    rm(pr.subset)

    print('Lon loop time')
    print(proc.time()-ltm)
  }##Longitude Loop
  nc_close(tasmax.nc)
  nc_close(tasmin.nc)
  nc_close(pr.nc)

  print('Removing lat band files')
  file.remove(paste0(tmp.dir,"/",tasmax.input))
  file.remove(paste0(tmp.dir,"/",tasmin.input))
  file.remove(paste0(tmp.dir,"/",pr.input))

}##Latitude File Loop

##Move back
write.dir <- paste0(base.dir,gcm)
move.back <- paste("rsync -av ",tmp.dir,scenario,"/",type," ",write.dir,"/",scenario ,sep='')
##print(move.back)
##system(move.back)

file.copy(from=paste0(tmp.dir,scenario,"/",type,"/"),to=paste0(write.dir,'/',scenario,'/'),overwrite=TRUE,recursive=TRUE)


if (type=='degree_days') {
  for (d in seq_along(degree.names)) {
    nc_close(dd.ncs[[d]]) 
  }
}

if (type=='annual') {
  for (d in seq_along(ann.names)) {
    nc_close(ann.ncs[[d]])
  }
}

if (type=='seasonal') {
  for (d in seq_along(seas.names)) {
    nc_close(seas.ncs[[d]])
  }
}

if (type=='monthly') {
  for (d in seq_along(mon.names)) {
    nc_close(mon.ncs[[d]])
  }
}

if (type=='annual_extremes') {
  for (d in seq_along(ext.names)) {
    nc_close(ext.ncs[[d]])
  }
}

if (type=='climdex') {
  for (d in seq_along(clim.names)) {
    nc_close(clim.ncs[[d]])
  }
}

##Rprof(NULL)


clean.up <- paste("rm ",tmp.dir,scenario,"/",type,"/*",sep='')
print(clean.up)
system(clean.up)

clean.up <- paste("rmdir ",tmp.dir,scenario,"/",type,sep='')
print(clean.up)
system(clean.up)

print("Total Elapsed Time:")
print(proc.time()-ptm)

