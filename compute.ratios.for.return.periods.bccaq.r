##Script to calculate return periods from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

library(ncdf4)
library(PCICt)
library(extRemes)
library(ismev)
library(udunits2)
library(zoo)

get.annual.data <- function(data,yearly.fac,var.name) {

  fx <- switch(var.name,
               tasmax=max,
               tasmin=min,
               pr=max)
  
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

  return(ts.to.fit)
  
}

calc.rp.ratios <- function(past.data,proj.data,
                           past.yearly.fac,proj.yearly.fac,
                           var.name,rperiod) {
  
  if (sum(is.na(past.data)) == length(past.data)) {
    return(NA)
  } else {
    past.to.fit <- get.annual.data(past.data,past.yearly.fac,var.name)
    proj.to.fit <- get.annual.data(proj.data,proj.yearly.fac,var.name)
    
    past.ts.fit <- fevd(past.to.fit,type='GEV')
    proj.ts.fit <- fevd(proj.to.fit,type='GEV')

    past.ts.old <- gev.fit(past.to.fit,show=FALSE)
    proj.ts.old <- gev.fit(proj.to.fit,show=FALSE)

    past.ts.fit$results$par <- past.ts.old$mle
    names(past.ts.fit$results$par) <- c('location','scale','shape')

    proj.ts.fit$results$par <- proj.ts.old$mle
    names(proj.ts.fit$results$par) <- c('location','scale','shape')

    ##past.rps <- return.level(past.ts.fit,return.period=as.numeric(rperiod),make.plot=F)
    past.rps <- return.level(past.ts.fit,return.period=as.numeric(rperiod),make.plot=F)

    new.rp <- 1/(1-pextRemes(proj.ts.fit,as.numeric(past.rps)))
    ##rp.ratio <- round(as.numeric(rperiod)/new.rp,1)

    ##if (rp.ratio > 10) {
    ##   browser()
    ##}
    if (sum(new.rp > 100) > 0) {
       browser()
    }
    new.rp[new.rp>100] <- 100
    new.rp[is.infinite(new.rp)] <- NA
    new.rp[new.rp<0.2] <- 0.2
    ##return(rp.ratio)
    return(new.rp)
  }
}


make.new.netcdf.file <- function(gcm,rcm=NULL,scenario,var.name,rperiod,
                                 past.file,write.file,
                                 data.dir,write.dir) {

  rp.name <- paste('rp.',rperiod,sep='')

  ##--------------------------------------------------------------
  nc <- nc_open(paste0(data.dir,'/',past.file),write=FALSE)
  
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
  ncatt_put(hist.nc,varid=rp.name,attname='units',attval='ratio')
  
  nc_close(hist.nc)  
  nc_close(nc)  

}

get.time <- function(nc) {

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time') 
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + ncvar_get(nc,'time')*86400
  return(var.dates)
}

date.bounds <- function(var.dates,interval) {
  bnds <- strsplit(interval,'-')[[1]]
  yst <- head(grep(bnds[1],var.dates),1)
  yen <- tail(grep(bnds[2],var.dates),1)
  if (length(yen)==0) {
     yen <- length(var.dates)
  }
  return(c(yst,yen))
}



rp.for.model <- function(gcm,rcm=NULL,scenario,var.name,rperiod,
                         past.file,proj.file,write.file,
                         past.int,proj.int,
                         data.dir,write.dir) {
                         

  rp.name <- paste('rp.',rperiod,sep='')  

  nc <- nc_open(paste(write.dir,write.file,sep=''),write=TRUE)
  hist.nc <- nc_open(paste0(data.dir,'/',past.file),write=FALSE)
  proj.nc <- nc_open(paste0(data.dir,'/',proj.file),write=FALSE)
  var.units <- ncatt_get(hist.nc,var.name,'units')$value
  print(var.units)    
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat) 

  past.dates <- get.time(hist.nc)
  proj.dates <- get.time(proj.nc)

  past.bnds <- date.bounds(past.dates,past.int)
  pst <- past.bnds[1]
  pen <- past.bnds[2]
  
  proj.bnds <- date.bounds(proj.dates,proj.int)
  fst <- proj.bnds[1]
  fen <- proj.bnds[2]

  
  past.yearly.fac <- as.factor(format(past.dates,'%Y'))[pst:pen]
  proj.yearly.fac <- as.factor(format(proj.dates,'%Y'))[fst:fen]

  print(n.lon)
  for (i in 247:n.lon) {
#    print(paste('i = ',i,sep=''))
      print(paste(gcm,'-',rcm,' i= ',i,sep=''))
      past.subset <- ncvar_get(hist.nc,var.name,start=c(i,1,pst),count=c(1,-1,(pen-pst+1)))
      proj.subset <- ncvar_get(proj.nc,var.name,start=c(i,1,fst),count=c(1,-1,(fen-fst+1)))

      past.list <- vector(mode='list',length=n.lat)
      proj.list <- vector(mode='list',length=n.lat)
      for (j in 1:n.lat) {
        past.list[[j]] <- past.subset[j,]
        proj.list[[j]] <- proj.subset[j,]
      }

      rp.values <- mapply(FUN=calc.rp.ratios,past.list,proj.list,MoreArgs=list(past.yearly.fac,proj.yearly.fac,var.name=var.name,rperiod))  

##      ncvar_put(nc,varid=rp.name,vals=rp.values,
##                   start=c(i,1,1),count=c(1,n.lat,1))            
  }
  
  nc_close(hist.nc)
  nc_close(nc)
}

##**************************************************************************************

##

rcp26.list <- c('CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-ES',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')
gcm.list <- c('CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G',
              'HadGEM2-CC',
              'HadGEM2-ES',
              'inmcm4',
              'MIROC5',
              'MPI-ESM-LR',
              'MRI-CGCM3')
##

gcm.list <- c('ACCESS1-0',
              'CanESM2',
              'CCSM4',
              'CNRM-CM5',
              'CSIRO-Mk3-6-0',
              'GFDL-ESM2G')
##              'HadGEM2-CC',
##              'HadGEM2-ES',
##              'inmcm4',
##              'MIROC5',
##              'MPI-ESM-LR',
##              'MRI-CGCM3')

gcm.list <- 'CCSM4'

###--------------------------------------------------------------------

run.bccaq.gcms.rp <- function() {
 
  var.name <- 'pr'
  scenario <- 'rcp85'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  rperiod <- '20'
  region <- 'van_whistler'

  data.dir <- paste('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/',sep='') 
  rp.dir <- paste('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/',scenario,'/return_periods/',sep='') 


  tmp.base <- paste('/local_temp/ssobie/',region,'/',var.name,'/',sep='')
  tmp.rp <- paste(tmp.base,scenario,'/return_periods/',sep='')
  
  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    rcm <- NULL
    write.dir <- paste0(rp.dir,gcm)

    tmp.dir <- paste(tmp.base,gcm,sep='')
    if (!file.exists(tmp.dir))
       dir.create(tmp.dir,recursive=TRUE)

    write.rp <- paste(tmp.rp,gcm,'/',sep='')
    if (!file.exists(write.rp))
      dir.create(write.rp,recursive=TRUE)

    var.files <- list.files(path=paste(data.dir,gcm,'/',sep=''),pattern=paste(var.name,'_gcm',sep=''))

    ##-------------------------------------------------    
    var.past.file <- var.files[grep('1951-2000',var.files)] ##var.files[grep(past.int,var.files)]
    var.proj.file <- var.files[grep('2001-2100',var.files)] ##var.files[grep(proj.int,var.files)]

    file.split <- strsplit(var.past.file,'_')[[1]]
    run <- file.split[grep('r*i1p1',file.split)]

    write.name <- paste(var.name,'_RP',rperiod,'_NEW_INTERVAL_BCCAQ_PRISM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')

    move.to <- paste("rsync -av ",data.dir,gcm,"/",var.past.file," ",tmp.dir,sep='')
    system(move.to)
    move.to <- paste("rsync -av ",data.dir,gcm,"/",var.proj.file," ",tmp.dir,sep='')
    system(move.to)


    make.new.netcdf.file(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,write.name,
                         tmp.dir,write.rp) ##data.dir,write.dir)

    print('made new file')
    test <- rp.for.model(gcm,rcm,scenario,var.name,rperiod,
                         var.past.file,var.proj.file,write.name,
                         past.int,proj.int,
                         tmp.dir,write.rp)

    move.back <- paste("rsync -av ",write.rp," ",write.dir,sep='')
    print(move.back)
    system(move.back)

    clean.up <- paste("rm ",tmp.dir,"/",var.past.file,sep='')
    system(clean.up)
    clean.up <- paste("rm ",tmp.dir,"/",var.proj.file,sep='')
    system(clean.up)

  }
}

##-------------------------------------------------------------------------

run.bccaq.gcms.rp()

