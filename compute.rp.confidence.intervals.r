##Script to calculate return periods from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file
##Modified to add the confidence intervals to the resulting netcdf file

library(ncdf4)
library(PCICt)
library(extRemes)
library(ismev)
library(udunits2)
library(zoo)

calc.return.periods <- function(data,yearly.fac,var.name,rperiod) {
  
  if (sum(is.na(data)) == length(data)) {
    return(list(ci=c(NA,NA,NA),
                rp=NA))
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
    ts.rp <- return.level(ts.fit,return.period=as.numeric(rperiod),make.plot=F)
    ts.rps <- c(NA,as.numeric(ts.rp),NA)

    tryCatch({
      ##print('In try')
      ts.rps <- ci(ts.fit,alpha=0.05,type=c('return.level','parameter'),return.period=as.numeric(rperiod))
    }, warning=function(war) {
      message('Warning from confidence intervals')      
    }, error=function(err) {
      print('Error in CI')
    }, finally={
      ##print(ts.rps)
    })

    ##This gives a 3 element vector with lower bound, rp value and upper bound
    rv <- as.numeric(ts.rps)
    rv <- list(ci=as.numeric(ts.rps),
               rp=as.numeric(ts.rp))
    
    if (var.name=='tasmin') {
      ##rv <- -as.numeric(ts.rps) ##$return.level
      rv <- list(ci=-as.numeric(ts.rps)[c(3,2,1)],
                 rp=-as.numeric(ts.rp))
    }
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
  t.geog <- ncdim_def('time', 'RP', c(1,2,3),
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
  ncatt_put(hist.nc,varid=rp.name,attname='units',attval=var.units)   ##'kg m-2 d-1')
  
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
  bnds <- strsplit(interval,'-')[[1]]
  
  st <- head(grep(bnds[1],var.dates),1)
  en <- tail(grep(bnds[2],var.dates),1)
  #st <- 1
  #en <- length(var.dates)

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
      rp.ci <- lapply(rp.values,function(x) {return(x$ci)})
      rp.rp <- lapply(rp.values,function(x) {return(x$rp)})
      print('rp.values')
      
      ##CI
      print(max(unlist(rp.ci),na.rm=T))
      rp.ci.checked <- lapply(rp.ci,check.rp.outliers,var.name)
      rp.ci.write <- matrix(unlist(rp.ci.checked),nrow=n.lat,ncol=3,byrow=TRUE)

      ##RP
      print(max(unlist(rp.rp),na.rm=T))
      rp.rp.checked <- lapply(rp.rp,check.rp.outliers,var.name)
      rp.rp.write <- matrix(unlist(rp.rp.checked),nrow=n.lat,ncol=1,byrow=TRUE)
      ncvar_put(nc,varid=rp.name,vals=rp.ci.write,
                   start=c(i,1,1),count=c(1,n.lat,3))            
  }
  
  nc_close(hist.nc)
  nc_close(nc)
}

##**************************************************************************************




###--------------------------------------------------------------------

run.bccaq.gcms.rp <- function() {

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
 
  var.name <- 'tasmin'
  scenario <- 'rcp85'
  past.int <- '1971-2000'
  proj.int <- '2071-2100'
  rperiod <- '20'

  data.dir <- paste('/storage/data/scratch/ssobie/bccaq_gcm_bc_subset/',sep='') 
  rp.dir <- paste('/storage/data/scratch/ssobie/bccaq_gcm_bc_subset/',scenario,'/return_periods/',sep='')
  
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    write.dir <- paste(rp.dir,gcm,'/',sep='')
    
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    var.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern=paste(var.name,'_day',sep=''),full.name=TRUE)
    scen.files <- var.files[grep(scenario,var.files)]    
    ##-------------------------------------------------    
    var.past.file <- scen.files[grep('1951-2000',scen.files)]

    file.split <- strsplit(var.past.file,'_')[[1]]
    run <- file.split[grep('r*i1p1',file.split)]

    write.hist.name <- paste(var.name,'_RPCI',rperiod,'_BCCAQ_GCM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')

    if (1==0) {
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,write.hist.name,
                         data.dir,write.dir,interval=past.int,canada=TRUE)
    } 
    if (1==1) {
    ##-------------------------------------------------
    var.proj.file <- scen.files[grep('2001-2100',scen.files)]
    write.proj.name <- paste(var.name,'_RPCI',rperiod,'_BCCAQ_GCM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')

    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         var.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         var.proj.file,write.proj.name,
                         data.dir,write.dir,interval=proj.int,canada=TRUE)
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
    data.dir <- '/storage/data/climate/downscale/CMIP5/BCCAQ/'
    rp.dir <- '/storage/data/climate/downscale/CMIP5/BCCAQ/return_periods/'
  }
  if (country=='America') {
    data.dir <- '/datasets/climate-downscale-CMIP5/nobackup/tmp_from_westgrid/'
    rp.dir <- '/storage/data/climate/downscale/CMIP5/BCCAQ/return_periods/america/'
  }
    
  gcm.list <- 'MPI-ESM-LR'
  
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    write.dir <- paste(rp.dir,gcm,'/',sep='')
    
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
    gcm.files <- list.files(path=data.dir,pattern=paste('BCCAQ\\+ANUSPLIN300\\+',model,'_',sep=''),full.name=TRUE)
    var.file <- gcm.files[grep(scenario,gcm.files)]

    ##-------------------------------------------------    
    file.split <- strsplit(var.file,'_')[[1]]
    run <- file.split[grep('r*i1p1',file.split)]

    if (1==1) {
      print('historical')
      write.hist.name <- paste(var.name,'_RP',rperiod,'_BCCAQ_GCM_',country,'_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')
      
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
      write.proj.name <- paste(var.name,'_RP',rperiod,'_BCCAQ_GCM_',country,'_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
      
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

  data.dir <- '/storage/data/scratch/ssobie/bccaq_gcm/rcp45/gcm/' 
  rp.dir <- '/storage/data/scratch/ssobie/bccaq_gcm/rcp45/gcm/return_periods_from_gcm_scale/' 

  var.name <- 'pr'
  scenario <- 'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2071-2100'
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
    write.hist.name <- paste('pr_RPCI',rperiod,'_GCM_',gcm,'_',scenario,'_gcm_scale_',past.int,'.nc',sep='')
    
    ##make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
    ##                     pr.past.file,write.hist.name,
    ##                     data.dir,write.dir)
    ##print('made new file')
    ##test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
    ##                     pr.past.file,write.hist.name,
    ##                     data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    ##write.proj.name <- paste('pr_RP',rperiod,'_',gcm,'_',scenario,'_',run,'_gcm_scale_',proj.int,'.nc',sep='')
    write.proj.name <- paste('pr_RPCI',rperiod,'_GCM_',gcm,'_',scenario,'_gcm_scale_',proj.int,'.nc',sep='')
   
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
  scenario <- 'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '10'

  data.dir <- paste('/storage/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/',sep='')
  ##rp.dir <- paste('/storage/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/return_periods/',sep='')
  rp.dir <- paste('/storage/data/scratch/ssobie/gcm/return_periods/',sep='')

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

    run <- strsplit(pr.past.file,'_')[[1]][6]
    
    write.hist.name <- paste('pr_RPCI',rperiod,'_GCM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')

    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    write.proj.name <- paste('pr_RPCI',rperiod,'_GCM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.proj.file,write.proj.name,
                         data.dir,write.dir)

  }
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
                   c('HADCM3','HRM3'),
                   c('HADCM3','MM5I'))
  
  data.dir <-  '/storage/data/scratch/ssobie/bccaq_rcm_bc_subset/' ##'/storage/data/projects/rci/data/stat.downscaling/scaling_comparison/bccaq_gcm/gcm/'
  rp.dir <-  '/storage/data/scratch/ssobie/bccaq_rcm/return_periods/'  ##'/storage/data/projects/rci/data/stat.downscaling/scaling_comparison/bccaq_gcm/gcm/return_periods/'

  var.name <- 'pr'
  scenario <- 'sresa2'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '20'
  
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
    write.hist.name <- paste('pr_RPCI',rperiod,'_BCCAQ_RCM_',gcm,'_',rcm,'_',scenario,'_',past.int,'.nc',sep='')
   
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    write.proj.name <- paste('pr_RPCI',rperiod,'_BCCAQ_RCM_',gcm,'_',rcm,'_',scenario,'_',proj.int,'.nc',sep='')
   
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
  
  data.dir <-  paste('/storage/data/scratch/ssobie/bccaq_rcm/',scale,'/',sep='') 
  rp.dir <-  paste('/storage/data/scratch/ssobie/bccaq_rcm/',scale,'/return_periods_from_',scale,'_scale/',sep='')  

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

  ##data.dir <- paste('/storage/data/projects/rci/data/stat.downscaling/scaling_comparison/rcm/gcm/',sep='')
  data.dir <- paste('/storage/data/scratch/ssobie/rcm/',sep='')
  rp.dir <- paste('/storage/data/scratch/ssobie/rcm/return_periods/',sep='')
  
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
    write.hist.name <- paste('pr_RPCI',rperiod,'_RCM_',gcm,'_',rcm,'_sresa2_rcm_scale_',past.int,'.nc',sep='')

    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         pr.past.file,write.hist.name,
                         data.dir,write.dir)

    ##-------------------------------------------------
    pr.proj.file <- pr.files[grep(proj.int,pr.files)]
    write.proj.name <- paste('pr_RPCI',rperiod,'_RCM_',gcm,'_',rcm,'_sresa2_rcm_scale_',proj.int,'.nc',sep='')
   
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

  scale <- 'raw'
  
  data.dir <-  paste('/storage/data/scratch/ssobie/anusplin/',sep='') 
  rp.dir <-  paste('/storage/data/scratch/ssobie/anusplin/return_periods/',sep='')  

  var.name <- 'pr'
  scenario <- 'v2013'
  past.int <- '1951-2000'
  rperiod <- '20'
 
  gcm <- 'ANUSPLIN'
  rcm <- NULL
  write.dir <- rp.dir
  
  if (!file.exists(write.dir))
    dir.create(write.dir,recursive=TRUE)
  
  pr.past.file <- paste(data.dir,'pr_day_ANUSPLIN_observation_v2013_1951-2000.nc',sep='')
  ##pr.past.file <- paste(data.dir,'pr_day_ANUSPLIN_observation_v2013_rcm_scale_1951-2000.nc',sep='')
  write.hist.name <- paste('pr_RPCI',rperiod,'_ANUSPLIN_',scenario,'_',scale,'_scale_',past.int,'.nc',sep='')

  make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir)
  print('made new file')
  test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                       pr.past.file,write.hist.name,
                       data.dir,write.dir,interval=past.int,canada=TRUE)

}


##-------------------------------------------------------------------------

run.bccaq.prism.rp <- function() {

##'ACCESS1-0',
##              'CanESM2',
##              'CCSM4',
##              'CNRM-CM5',
##              'CSIRO-Mk3-6-0',
##              'GFDL-ESM2G',
##              'HadGEM2-CC',
##              'HadGEM2-ES',

gcm.list <- c('inmcm4',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')

  var.name <- 'tasmin'
  scenario <- 'rcp85'
  past.int <- '1971-2000'
  proj.int <- '2071-2100'
  rperiod <- '10'
  region <- 'van_whistler'

  data.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/')
  write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/',scenario,'/return_periods/')
  tmp.base <- paste('/local_temp/ssobie/',region,'/',var.name,'/',sep='')

  tmp.rp <- paste(tmp.base,scenario,'/return_periods/',sep='')
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    tmp.dir <- paste(tmp.base,gcm,sep='')
    if (!file.exists(tmp.dir))
       dir.create(tmp.dir,recursive=TRUE)
    move.to <- paste("rsync -av ",data.dir,gcm,"/",var.name,"_gcm_prism_BCCAQ_",gcm,"* ",tmp.dir,sep='')
    print(move.to)
    system(move.to)

    write.rp <- paste(tmp.rp,gcm,'/',sep='')    
    if (!file.exists(write.rp))
      dir.create(write.rp,recursive=TRUE)
    
    var.files <- list.files(path=paste(tmp.dir,'/',sep=''),pattern=paste(var.name,'_gcm_prism',sep=''),full.name=TRUE)

    ##-------------------------------------------------    
    ##Past File                                                      
    var.past.file <- var.files[grep('1951-2000',var.files)]
    file.split <- strsplit(var.past.file,'_')[[1]]
    run <- file.split[grep('r*i1p1',file.split)]
    write.hist.name <- paste(var.name,'_RPCI',rperiod,'_BCCAQ_PRISM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')

    if (1==0) {
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,write.hist.name,
                         tmp.dir,write.rp)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,write.hist.name,
                         tmp.dir,write.rp,interval=past.int,canada=TRUE)
    } 
    if (1==1) {
    ##-------------------------------------------------
    ##Future File
    var.proj.file <- var.files[grep('2001-2100',var.files)]
    write.proj.name <- paste(var.name,'_RPCI',rperiod,'_BCCAQ_PRISM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         var.proj.file,write.proj.name,
                         tmp.dir,write.rp)
    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         var.proj.file,write.proj.name,
                         tmp.dir,write.rp,interval=proj.int,canada=TRUE)
  }

  move.back <- paste("rsync -av ",tmp.rp,gcm," ",write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up <- paste("rm ",tmp.dir,"/",var.name,"_gcm_prism_BCCAQ_",gcm,"* " ,sep='')
  print(clean.up)
  system(clean.up)

  clean.up.dd <- paste("rm ",write.rp,var.name,"_RPCI*" ,sep='')
  print(clean.up.dd)
  system(clean.up.dd)


  }##GCM loop
}

##-------------------------------------------------------------------------


run.bccaq.prism.rp() 
