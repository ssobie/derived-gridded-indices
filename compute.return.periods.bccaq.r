##Script to calculate return periods from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file
library(zoo)
library(ncdf4)
library(PCICt)
library(extRemes)
library(ismev)
library(udunits2)


calc.return.periods <- function(data,yearly.fac,var.name,rperiod) {
  
  if (sum(is.na(data)) == length(data)) {
    return(NA)
  } else {
    fx <- switch(var.name,
                 tasmax=max,
                 tasmin=min,
                 pr=max)
    
##    ts.sum <- tapply(data,list(yearly.fac),function(x) {return(rollapply(x,width=3,FUN=sum))})
##    ts.yearly <- unlist(lapply(ts.sum,max,na.rm=T)) 
    
    ts.yearly <- tapply(data,list(yearly.fac),fx,na.rm=T)
    inf.flag <- is.infinite(ts.yearly)
    ts.yearly[inf.flag] <- NA

    na.flag <- is.na(ts.yearly)
    if(sum(na.flag)>0) {
      ts.to.fit <- as.vector(ts.yearly[-which(na.flag)])
    } else {
      ts.to.fit <- as.vector(ts.yearly)
    }
    if (var.name=='tasmin') {
      ts.to.fit <- -ts.to.fit
      u.len <- length(unique(ts.to.fit))
      f.len <- length(ts.to.fit)
      if (u.len < (f.len*2/3))
        ts.to.fit <- jitter(ts.to.fit,amount=3)
    }      

    ts.fit <- fevd(ts.to.fit,type='GEV')
    ts.old <- gev.fit(ts.to.fit,show=FALSE)
    
    ts.fit$results$par <- ts.old$mle
    names(ts.fit$results$par) <- c('location','scale','shape')

    ts.rps <- return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)

    rv <- as.numeric(ts.rps)

    if (var.name=='tasmin')
      rv <- -as.numeric(ts.rps) ##$return.level
    return(rv)
  }
}


make.new.netcdf.file <- function(gcm,rcm=NULL,scenario,var.name,rperiod,
                                 var.file,write.file,
                                 data.dir,write.dir) {

  rp.name <- paste('rp.',rperiod,sep='')

  ##--------------------------------------------------------------
  nc <- nc_open(var.file,write=FALSE)
  
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time') 
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400
  years <- unique(format(time.series,'%Y'))

  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  
  
  lon.atts <- ncatt_get(nc,'lon')
  lat.atts <- ncatt_get(nc,'lat')
  global.atts <- ncatt_get(nc,varid=0)
  
  var.atts <- ncatt_get(nc,var.name)
  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##--------------------------------------------------------------
  ##Create new netcdf file

  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', 'RP', 1,
                      unlim=TRUE, calendar=time.calendar)
  
  var.geog <- ncvar_def(rp.name, units=var.atts$units, dim=list(x.geog, y.geog, t.geog),
                        missval=var.atts[['_FillValue']])
  
  hist.nc <- nc_create(paste(write.dir,write.file,sep=''), var.geog)
  
  ##Loop over subsets of the time series
  ##Past file first
  global.names <- names(global.atts)
  for (g in 1:length(global.atts)) 
    ncatt_put(hist.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ##Time attributes
  ncatt_put(hist.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(hist.nc,varid='time',attname='long_name',attval='Years')
  ncatt_put(hist.nc,varid='time',attname='standard_name',attval='Years')
  ncatt_put(hist.nc,varid='time',attname='calendar',attval=time.calendar)  
  
  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))  
    ncatt_put(hist.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])
  
  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))  
    ncatt_put(hist.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])

  ##Climdex Attributes
  var.names <- names(var.atts)
  for (j in 1:length(var.atts))
    ncatt_put(hist.nc,varid=rp.name,attname=var.names[j],attval=var.atts[[j]])

  var.units <- ncatt_get(nc,var.name,'units')$value
  ncatt_put(hist.nc,varid=rp.name,attname='units',attval='kg m-2 d-1')
  
  nc_close(hist.nc)  
  nc_close(nc)  

}

check.model.outliers <- function(data,var.name) {
  rv <- data
  if (var.name=='tasmax')
    flags <- which(data > 75)
  if (var.name=='tasmin')
    flags <- which(data < -75)
  if (var.name=='pr')
    flags <- which(data > 500)
  if (length(flags)!=0)
    rv[flags] <- NA

  return(rv)
}


check.rp.outliers <- function(data,var.name) {
  rv <- data
  if (var.name=='tasmax')
    flags <- which(data > 75)
  if (var.name=='tasmin')
    flags <- which(data < -75)
  if (var.name=='pr')
    flags <- which(data > 500)
  if (length(flags)!=0)
    rv[flags] <- 1111

  return(rv)
}


rp.for.model <- function(gcm,rcm=NULL,scenario,var.name,rperiod,
                         var.file,write.file,
                         data.dir,write.dir,
                         interval=NULL,canada=FALSE) {

  rp.name <- paste('rp.',rperiod,sep='')  

  nc <- nc_open(paste(write.dir,write.file,sep=''),write=TRUE)
  hist.nc <- nc_open(var.file,write=FALSE)
  var.units <- ncatt_get(hist.nc,var.name,'units')$value
  print(var.units)    
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat)
  
  time.atts <- ncatt_get(hist.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(hist.nc,'time') 
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  
  var.dates <- origin.pcict + ncvar_get(hist.nc,'time')*86400

  st <- 1
  en <- length(var.dates)

  if (canada) {
    yr.bnds <- strsplit(interval,'-')[[1]]
    st <- grep(yr.bnds[1],var.dates)[1]
    en <- tail(grep(yr.bnds[2],var.dates),1)    
    if (length(en) == 0)
      en <- length(var.dates)
    print(paste(yr.bnds[1],' at ',st,sep=''))
    print(paste(yr.bnds[2],' at ',en,sep=''))
  }

  yearly.fac <- as.factor(format(var.dates[st:en],'%Y'))

  print(n.lon)
  for (i in 1:n.lon) {
#    print(paste('i = ',i,sep=''))
      print(paste(gcm,'-',rcm,' i= ',i,sep=''))
      var.subset <- ncvar_get(hist.nc,var.name,start=c(i,1,st),count=c(1,-1,(en-st+1)))
      if (sum(is.na(var.subset)) != length(data))
        var.subset <- check.model.outliers(var.subset,var.name)

      var.list <- list()
      for (j in 1:n.lat) {
        var.list[[j]] <- var.subset[j,]
      }
      if (var.units == 'kg m-2 s-1')
        var.list <- lapply(var.list,ud.convert,var.units,'kg m-2 d-1')

      rp.values <- lapply(var.list,calc.return.periods,yearly.fac,var.name=var.name,rperiod)      
      print(max(unlist(rp.values),na.rm=T))
      rp.checked <- lapply(rp.values,check.rp.outliers,var.name)
      print(rp.checked)

      ncvar_put(nc,varid=rp.name,vals=rp.checked,
                   start=c(i,1,1),count=c(1,n.lat,1))            
  }
  
  nc_close(hist.nc)
  nc_close(nc)
}

##**************************************************************************************
gcm.list <- c('CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-ES',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')


gcm.list <- c('ACCESS1-0',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-CC',
              'HadGEM2-ES',
              'inmcm4',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')

gcm.list <- c('inmcm4',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')


###--------------------------------------------------------------------

run.bccaq.gcms.rp <- function() {
 
  var.name <- 'pr'
  scenario <- 'rcp85'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '10'

  data.dir <- paste('/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_bc_subset/',scenario,'/',sep='') 
  rp.dir <- paste('/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm/',scenario,'/return_periods/',sep='')
  
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    write.dir <- paste(rp.dir,gcm,'/',sep='')
    
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    var.files <- list.files(path=paste(data.dir,gcm,'/',sep=''),pattern=paste(var.name,'_day',sep=''),full.name=TRUE)

    ##-------------------------------------------------    
    var.past.file <- var.files[grep(past.int,var.files)]
    run <- strsplit(var.past.file,'_')[[1]][9]

    write.hist.name <- paste(var.name,'_RP',rperiod,'_BCCAQ_GCM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')

    if (1==1) {
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,write.hist.name,
                         data.dir,write.dir)
    } 
    if (1==0) {
    ##-------------------------------------------------
    var.proj.file <- var.files[grep(proj.int,var.files)]
    write.proj.name <- paste(var.name,'_RP',rperiod,'_BCCAQ_GCM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         var.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         var.proj.file,write.proj.name,
                         data.dir,write.dir)
  }
  }
}

run.country.bccaq.gcms.rp <- function() {

  ptm <- proc.time()

  country <- 'Canada'
  var.name <- 'tasmin'
  scenario <- 'rcp85'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '20'

  if (country=='Canada') {
    data.dir <- '/storage/data/climate/downscale/BCCAQ2/nobackup/CanESM2/tasmin/'
    rp.dir <- '/storage/data/projects/downscale-idf/BCCAQ2/return_periods/'
  }
  if (country=='America') {
    data.dir <- '/datasets/climate-downscale-CMIP5/nobackup/CONUS/'
    rp.dir <- '/home/data/climate/downscale/CMIP5/BCCAQ/return_periods/america/'
  }
    
  gcm.list <- 'CanESM2'
  
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    write.dir <- paste(rp.dir,gcm,'/',sep='')
    
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    ##gcm.files <- list.files(path=data.dir,pattern=paste('BCCAQ\\+ANUSPLIN300\\+',model,'_',sep=''),full.name=TRUE)
    ##gcm.files <- list.files(path=data.dir,pattern=paste('bccaq\\+',model,'_',sep=''),full.name=TRUE)
    
    ##var.file <- gcm.files[grep(scenario,gcm.files)]
    var.file <- paste(data.dir,'tasmin_output.nc',sep='')
    ##-------------------------------------------------    
    if (country=='Canada')
      run <- 'r1i1p1' ##strsplit(var.file,'_')[[1]][5]
    if (country=='America')
      run <- strsplit(var.file,'_')[[1]][7]

    if (1==1) {
      print('historical')
      write.hist.name <- paste(var.name,'_RP',rperiod,'_BCCAQ2_GCM_',country,'_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')
      
      make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                           var.file,write.hist.name,
                           data.dir,write.dir)
      print('made new file')
      test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                           var.file,write.hist.name,
                           data.dir,write.dir,
                           past.int,canada=TRUE)
    }
    ##-------------------------------------------------
    if (1==1) {
      print('projection')
      ##    var.proj.file <- var.files[grep(proj.int,var.files)]
      write.proj.name <- paste(var.name,'_RP',rperiod,'_BCCAQ2_GCM_',country,'_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
      
      make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                           var.file,write.proj.name,
                           data.dir,write.dir)
      print('made new file')
      test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                           var.file,write.proj.name,
                           data.dir,write.dir,
                           proj.int,canada=TRUE)

    }
    print('Elapsed Time')    
    print(proc.time()-ptm)
    
  }
}

run.bccaq.gcm.scale.rp <- function() {

  data.dir <- '/home/data/scratch/ssobie/bccaq_gcm/rcp45/gcm/' 
  rp.dir <- '/home/data/scratch/ssobie/bccaq_gcm/rcp45/gcm/return_periods_from_gcm_scale/' 

  var.name <- 'pr'
  scenario <- 'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '20'
  
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    write.dir <- paste(rp.dir,gcm,'/',sep='')
     
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    pr.files <- list.files(path=paste(data.dir,gcm,'/',sep=''),pattern='pr_day',full.name=TRUE)

    ##-------------------------------------------------    
    pr.past.file <- pr.files[grep(past.int,pr.files)]
    run <- strsplit(pr.past.file,'_')[[1]][9]
    
    ##write.hist.name <- paste('pr_RP',rperiod,'_',gcm,'_',scenario,'_',run,'_gcm_scale_',past.int,'.nc',sep='')
    write.hist.name <- paste('pr_RP',rperiod,'_GCM_',gcm,'_',scenario,'_gcm_scale_',past.int,'.nc',sep='')
    
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    ##write.proj.name <- paste('pr_RP',rperiod,'_',gcm,'_',scenario,'_',run,'_gcm_scale_',proj.int,'.nc',sep='')
    write.proj.name <- paste('pr_RP',rperiod,'_GCM_',gcm,'_',scenario,'_gcm_scale_',proj.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)

  }
}


run.gcms.rp <- function() {

  var.name <- 'pr'
  scenario <- 'sresa2' ##'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '10'

  data.dir <- '/home/data/scratch/ssobie/cmip3_bc_subset/20c3m/' ##paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/',sep='')
  ##rp.dir <- paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/return_periods/',sep='')
  rp.dir <- paste('/home/data/scratch/ssobie/gcm/return_periods/',sep='')
  gcm.list <- 'CCCMA_CGCM3_1'
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    write.dir <- paste(rp.dir,gcm,'/',sep='')
     
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    pr.files <- list.files(path=paste(data.dir,tolower(gcm),'/',sep=''),pattern='pr_day',full.name=TRUE)
    browser()
    ##-------------------------------------------------    
    pr.past.file <- pr.files[grep(past.int,pr.files)]

    run <- strsplit(pr.past.file,'_')[[1]][6]
    
    write.hist.name <- paste('pr_RP',rperiod,'_GCM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')

    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    data.dir <- '/home/data/scratch/ssobie/cmip3_bc_subset/sresa2/'
    pr.files <- list.files(path=paste(data.dir,tolower(gcm),'/',sep=''),pattern='pr_day',full.name=TRUE)
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    write.proj.name <- paste('pr_RP',rperiod,'_GCM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)

  }
}

run.miroc4h.rp <- function() {

  var.name <- 'pr'
  scenario <- 'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2006-2035'
  rperiod <- '10'

  data.dir <- paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm_bc_subset/',scenario,'/',sep='')
  ##rp.dir <- paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/return_periods/',sep='')
  rp.dir <- paste('/home/data/scratch/ssobie/gcm/return_periods/',sep='')

  model <- 'MIROC4h'
  print(model)
  gcm <- model
  rcm <- NULL
  write.dir <- paste(rp.dir,gcm,'/',sep='')
     
  if (!file.exists(write.dir))
    dir.create(write.dir,recursive=TRUE)
  
  pr.files <- list.files(path=paste(data.dir,gcm,'/',sep=''),pattern='pr_day',full.name=TRUE)
  
  ##-------------------------------------------------    
  pr.past.file <- pr.files[grep(past.int,pr.files)]
  browser()
  run <- strsplit(pr.past.file,'_')[[1]][6]
  
  write.hist.name <- paste('pr_RP',rperiod,'_GCM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')

  make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir)
  print('made new file')
  test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir)
  
  ##-------------------------------------------------
  pr.proj.file <- pr.files[grep(proj.int,pr.files)]
  write.proj.name <- paste('pr_RP',rperiod,'_GCM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
  
  make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                       pr.proj.file,write.proj.name,
                       data.dir,write.dir)
  print('made new file')
  test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                       pr.proj.file,write.proj.name,
                       data.dir,write.dir)
    
}


run.canrcm4.rp <- function() {

  var.name <- 'pr'
  scenario <- 'rcp85'
  past.int <- '1950-2006'
  proj.int <- '2006-2100'
  rperiod <- '10'

  data.dir <- '/datasets/climate-downscale-idf-ec/data/CanRCM4-NAM-44/CCCma-CanESM2/sub_daily/'
  rp.dir <- paste('/home/data/scratch/ssobie/rcm/return_periods/',sep='')

  model <- 'CanRCM4'
  print(model)
  gcm <- model
  rcm <- NULL
  write.dir <- paste(rp.dir,gcm,'/',sep='')
     
  if (!file.exists(write.dir))
    dir.create(write.dir,recursive=TRUE)
  
  ##-------------------------------------------------    
  pr.past.file <- paste(data.dir,'pr_NAM-44_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_24hr_1950010101-2006010100.nc',sep='')
  run <- 'r1i1p1'
  
  write.hist.name <- paste('pr_RP',rperiod,'_RCM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')
  make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir)
  print('made new file')
  test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir)
  
  ##-------------------------------------------------
  pr.proj.file <- paste(data.dir,'pr_NAM-44_CCCma-CanESM2_rcp85_r1i1p1_CCCma-CanRCM4_24hr_2006010101-2100010100.nc',sep='')
  write.proj.name <- paste('pr_RP',rperiod,'_RCM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
  
  make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                       pr.proj.file,write.proj.name,
                       data.dir,write.dir)
  print('made new file')
  test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                       pr.proj.file,write.proj.name,
                       data.dir,write.dir)
    
}



###************************************************************************************************
###************************************************************************************************
##RCM data

run.bccaq.rcms.rp <- function() {

  rcm.list <- list(c('CCSM','CRCM'),
                   c('CCSM','MM5I'),
                   c('CCSM','WRFG'),
                   c('CGCM3','CRCM'),
                   c('CGCM3','RCM3'),
                   c('CGCM3','WRFG'),
                   c('GFDL','ECP2'),
                   c('GFDL','HRM3'),
                   c('GFDL','RCM3'),
                   c('HADCM3','HRM3'))
  rcm.list <- list(c('HADCM3','MM5I'))
  
  data.dir <-  '/home/data/scratch/ssobie/bccaq_rcm_bc_subset/' ##'/home/data/projects/rci/data/stat.downscaling/scaling_comparison/bccaq_gcm/gcm/'
  rp.dir <-  '/home/data/scratch/ssobie/bccaq_rcm/return_periods/'  ##'/home/data/projects/rci/data/stat.downscaling/scaling_comparison/bccaq_gcm/gcm/return_periods/'

  var.name <- 'pr'
  scenario <- 'sresa2'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '10'
  
  for (model in rcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- model[2]
    write.dir <- paste(rp.dir,gcm,'_',rcm,'/',sep='')
     
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    pr.files <- list.files(path=paste(data.dir,gcm,'_',rcm,'/',sep=''),pattern='pr_day',full.name=TRUE)

    ##-------------------------------------------------    
    pr.past.file <- pr.files[grep(past.int,pr.files)]
    write.hist.name <- paste('pr_RP',rperiod,'_BCCAQ_RCM_',gcm,'_',rcm,'_',scenario,'_',past.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    write.proj.name <- paste('pr_RP',rperiod,'_BCCAQ_RCM_',gcm,'_',rcm,'_',scenario,'_',proj.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    
  }
}


run.bccaq.rcms.agg.scales.rp <- function() {

  rcm.list <- list(c('CCSM','CRCM'),
                   c('CCSM','MM5I'),
                   c('CCSM','WRFG'),
                   c('CGCM3','CRCM'),
                   c('CGCM3','RCM3'),
                   c('CGCM3','WRFG'),
                   c('GFDL','ECP2'),
                   c('GFDL','HRM3'),
                   c('GFDL','RCM3'),
                   c('HADCM3','HRM3'),
                   c('HADCM3','MM5I'))

  scale <- 'rcm'
  
  data.dir <-  paste('/home/data/scratch/ssobie/bccaq_rcm/',scale,'/',sep='') 
  rp.dir <-  paste('/home/data/scratch/ssobie/bccaq_rcm/',scale,'/return_periods_from_',scale,'_scale/',sep='')  

  var.name <- 'pr'
  scenario <- 'sresa2'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '10'
  
  for (model in rcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- model[2]
    write.dir <- paste(rp.dir,gcm,'_',rcm,'/',sep='')
     
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    pr.files <- list.files(path=paste(data.dir,gcm,'_',rcm,'/',sep=''),pattern='pr_day',full.name=TRUE)

    ##-------------------------------------------------    
    pr.past.file <- pr.files[grep(past.int,pr.files)]
    write.hist.name <- paste('pr_RP',rperiod,'_BCCAQ_RCM_',gcm,'_',rcm,'_',scenario,'_',scale,'_scale_',past.int,'.nc',sep='')

    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    write.proj.name <- paste('pr_RP',rperiod,'_BCCAQ_RCM_',gcm,'_',rcm,'_',scenario,'_',scale,'_scale_',proj.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    browser()
  }
}

run.rcms.rp <- function() {

  rcm.list <- list(c('CCSM','CRCM'),
                   c('CCSM','MM5I'),
                   c('CCSM','WRFG'),
                   c('CGCM3','CRCM'),
                   c('CGCM3','RCM3'),
                   c('CGCM3','WRFG'),
                   c('GFDL','ECP2'),
                   c('GFDL','HRM3'),
                   c('GFDL','RCM3'),
                   c('HADCM3','HRM3'),
                   c('HADCM3','MM5I'))  

  var.name <- 'pr'
  scenario <- 'sresa2'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '20'

  ##data.dir <- paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/rcm/gcm/',sep='')
  data.dir <- paste('/home/data/scratch/ssobie/rcm/',sep='')
  rp.dir <- paste('/home/data/scratch/ssobie/rcm/return_periods/',sep='')
  
  for (model in rcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- model[2]
    write.dir <- paste(rp.dir,gcm,'_',rcm,'/',sep='')
     
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    pr.files <- list.files(path=paste(data.dir,tolower(gcm),'.',tolower(rcm),'/',sep=''),pattern='pr_day',full.name=TRUE) ##For RCMs
    ##pr.files <- list.files(path=paste(data.dir,gcm,'_',rcm,'/',sep=''),pattern='pr_day',full.name=TRUE) ##For GCMs

    ##-------------------------------------------------    
    pr.past.file <- pr.files[grep(past.int,pr.files)]
    write.hist.name <- paste('pr_RP',rperiod,'_RCM_',gcm,'_',rcm,'_sresa2_rcm_scale_',past.int,'.nc',sep='')

    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    write.proj.name <- paste('pr_RP',rperiod,'_RCM_',gcm,'_',rcm,'_sresa2_rcm_scale_',proj.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)

  }
}


##-------------------------------------------------------------------------

run.anusplin.rp <- function() {

  scale <- 'gcm'
  
  data.dir <-  paste('/home/data/scratch/ssobie/anusplin/gcm/',sep='') 
  rp.dir <-  paste('/home/data/scratch/ssobie/anusplin/gcm/return_periods/',sep='')  

  var.name <- 'pr'
  scenario <- 'v2013'
  past.int <- '1971-2000'
  rperiod <- '10'
 
  gcm <- 'ANUSPLIN'
  rcm <- NULL
  write.dir <- rp.dir
  
  if (!file.exists(write.dir))
    dir.create(write.dir,recursive=TRUE)
  
  pr.past.file <- paste(data.dir,'pr_day_ANUSPLIN_observation_v2013_gcm_scale_1951-2000.nc',sep='')
  write.hist.name <- paste('test_pr_RP',rperiod,'_ANUSPLIN_',scenario,'_',scale,'_scale_',past.int,'.nc',sep='')

  make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir)
  print('made new file')
  test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir)


}


##-------------------------------------------------------------------------

##Rprof('rp.bccaq.gcms.profile.out')
#run.bccaq.gcms.rp()
run.country.bccaq.gcms.rp()
##run.bccaq.rcms.agg.scales.rp()
##run.rcms.rp()
##Rprof(NULL)

