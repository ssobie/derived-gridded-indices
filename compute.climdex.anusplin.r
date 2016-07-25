##Script to calculate the climdex indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

##Modified to run just anusplin

library(ncdf4)
library(PCICt)
library(climdex.pcic)

get.climdex.info <- function(climdex.name) {
  
  climdex.names <- list(climdex.fd=c('fdETCCDI','Ann','days'),
                        climdex.su=c('suETCCDI','Ann','days'),
                        climdex.id=c('idETCCDI','Ann','days'),
                        climdex.tr=c('trETCCDI','Ann','days'),
                        climdex.gsl=c('gslETCCDI','Ann','days'),
                        climdex.txx=c('txxETCCDI','Mon','degC'),
                        climdex.tnx=c('tnxETCCDI','Mon','degC'),
                        climdex.txn=c('txnETCCDI','Mon','degC'),
                        climdex.tnn=c('tnnETCCDI','Mon','degC'),
                        climdex.tn10p=c('tn10pETCCDI','Mon','days'),
                        climdex.tx10p=c('tx10pETCCDI','Mon','days'),
                        climdex.tn90p=c('tn90pETCCDI','Mon','days'),
                        climdex.tx90p=c('tx90pETCCDI','Mon','days'),
                        climdex.wsdi=c('wsdiETCCDI','Ann','days'),
                        climdex.csdi=c('csdiETCCDI','Ann','days'),
                        climdex.dtr=c('dtrETCCDI','Mon','degC'),
                        climdex.rx1day=c('rx1dayETCCDI','Mon','mm'),
                        climdex.rx5day=c('rx5dayETCCDI','Mon','mm'),
                        climdex.sdii=c('sdiiETCCDI','Ann','mm d-1'),
                        climdex.r10mm=c('r10mmETCCDI','Ann','days'),
                        climdex.r20mm=c('r20mmETCCDI','Ann','days'),
                        climdex.cdd=c('cddETCCDI','Ann','days'),
                        climdex.cwd=c('cwdETCCDI','Ann','days'),
                        climdex.r95ptot=c('r95pETCCDI','Ann','mm'),
                        climdex.r99ptot=c('r99pETCCDI','Ann','mm'),
                        climdex.prcptot=c('prcptotETCCDI','Ann','mm'),
                        climdex.rnnmm=c('rnnmmETCCDI','Ann','mm'))                 
  rv <- climdex.names[[climdex.name]]
  return(rv)
}


get.climdex.value <- function(climdex.var,climdex.object) {
  
  rv <- eval(parse(text=paste(climdex.var,'(climdex.object)',sep='')))
  return(rv)
}


create.climdex.base.files <- function(climdex.name,gcm,rcm=NULL,type=NULL,
                                      past.int,
                                      data.dir,write.dir) {

  climdex.info <- get.climdex.info(climdex.name)
  climdex.var <- climdex.info[1]
  climdex.calendar <- climdex.info[2]
  climdex.units <- climdex.info[3]

  
  past.file <- list.files(path=data.dir,pattern='pr_gcm_prism',full.name=TRUE)
  if (!file.exists(write.dir))
    dir.create(write.dir,recursive=TRUE)
  
  ##write.clim.name <- paste(climdex.var,'_',tolower(climdex.calendar),'_anusplin_',past.int,'.nc',sep='')
  ##write.clim.name <- paste('r1mm_',tolower(climdex.calendar),'_anusplin_',past.int,'.nc',sep='')
  write.clim.name <- paste(climdex.var,'_',tolower(climdex.calendar),'_BCCAQ-PRISM_',gcm,'_rcp85_base_',past.int,'.nc',sep='')
  write.dir <- paste(write.dir,gcm,'/',sep='')

  nc <- nc_open(past.file,write=FALSE)
  
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)

  past.values <- ncvar_get(nc,'time')

  past.series <- format(past.origin + past.values*86400,'%Y-%m-%d')
  time.series <- past.series
  time.values <- as.Date(time.series) - time.start

  ##years <- seq(from=as.Date('1950-01-03'),by='year',to=as.Date('2005-12-03')) - time.start
  ##months <- seq(from=as.Date('1950-01-03'),by='month',to=as.Date('2005-12-03')) - time.start  

  years.ix <- grep('*-01-01',time.series)
  years <- time.values[years.ix]
  months.ix <- grep('[0-9]{4}-[0-9]{2}-01',time.series) ##grep('*-*-01',full.series)
  months <- time.values[months.ix] 
  
  dates <- switch(climdex.calendar,
                  Ann=years,
                  Mon=months)
  dates <- as.numeric(dates)

    ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  
  
  lon.atts <- ncatt_get(nc,'lon')
  lat.atts <- ncatt_get(nc,'lat')
  global.atts <- ncatt_get(nc,varid=0)
  
  pr.atts <- ncatt_get(nc,'pr')
  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, dates,
                      unlim=TRUE, calendar=time.calendar)

  var.geog <- ncvar_def(climdex.var, units=climdex.units, dim=list(x.geog, y.geog, t.geog),
                        missval=pr.atts[['_FillValue']])

  file.nc <- nc_create(paste(write.dir,write.clim.name,sep=''), var.geog)
  
  ##Loop over subsets of the time series
  ##Past file first
  global.names <- names(global.atts)
  for (g in 1:length(global.atts)) 
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ##Time attributes
  ncatt_put(file.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=time.calendar)  
  
  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))  
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])
  
  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))  
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])

  ##Climdex Attributes
  ncatt_put(file.nc,varid=climdex.var,attname='units',attval=climdex.units)
  ncatt_put(file.nc,varid=climdex.var,attname='_FillValue',attval=pr.atts[['_FillValue']])
  ncatt_put(file.nc,varid=climdex.var,attname='standard_name',attval=climdex.var)
  ncatt_put(file.nc,varid=climdex.var,attname='long_name',attval=climdex.var)

  nc_close(file.nc)
    
}


##---------------------------------------------------------------

climdex.for.model <- function(gcm,rcm=NULL,interval,type=NULL,
                              climdex.names,
                              past.int,
                              data.dir,write.dir) {

  pr.file <- list.files(path=data.dir,pattern='pr_gcm_prism',full.name=TRUE)
  tasmax.file <- list.files(path=data.dir,pattern='tasmax_gcm_prism',full.name=TRUE)
  tasmin.file <- list.files(path=data.dir,pattern='tasmin_gcm_prism',full.name=TRUE)

  hist.dir <- paste(write.dir,gcm,sep='')

  clim.files <- list.files(path=hist.dir,pattern='ETCCDI',full.name=TRUE)

  clim.ncs <- lapply(clim.files,nc_open,write=TRUE)

  ##--------------------------------------------------------------
  climdex.vars <- sort(unlist(lapply(lapply(climdex.names,get.climdex.info),function(x){return(x[1])})))
  
  pr.past.nc <- nc_open(pr.file,write=FALSE)
  tasmax.past.nc <- nc_open(tasmax.file,write=FALSE)
  tasmin.past.nc <- nc_open(tasmin.file,write=FALSE)

  pr.scale <- 1
  temp.offset <- 0
  pr.units <- ncatt_get(pr.past.nc,'pr')$units
  if (pr.units == 'kg m-2 s-1')
    pr.scale <- 86400
  tx.units <- ncatt_get(tasmax.past.nc,'tasmax')$units
  if (tx.units == 'K')
    temp.scale <- 273
  
  ##Attributes to retain
  lon <- ncvar_get(pr.past.nc,'lon')
  lat <- ncvar_get(pr.past.nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##Combine the dates
  time.atts <- ncatt_get(pr.past.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
 
  pr.past.values <- ncvar_get(pr.past.nc,'time')
  pr.dates <- past.origin + pr.past.values*86400

  tasmax.past.values <- ncvar_get(tasmax.past.nc,'time')
  tasmax.dates <- past.origin + tasmax.past.values*86400
 
  tasmin.past.values <- ncvar_get(tasmin.past.nc,'time')
  tasmin.dates <- past.origin + tasmin.past.values*86400

  ##--------------------------------------------------------------
  ##Compute climdex values and load into newly created climdex netcdf
  for (i in 1:n.lon) {
    print(paste('Lon: ',i,' in ',n.lon,sep=''))

    pr.subset <- ncvar_get(pr.past.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))*pr.scale
    tasmax.subset <- ncvar_get(tasmax.past.nc,'tasmax',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
    tasmin.subset <- ncvar_get(tasmin.past.nc,'tasmin',start=c(i,1,1),count=c(1,-1,-1))-temp.offset

    pr.list <- list()
    tasmax.list <- list()
    tasmin.list <- list()

      for (j in 1:n.lat) {
        tasmax.list[[j]] <- tasmax.subset[j,]
        tasmin.list[[j]] <- tasmin.subset[j,]
        pr.list[[j]] <- pr.subset[j,]
      }

      climdex.objects <- mapply(climdexInput.raw,tasmax.list,tasmin.list,pr.list,
                                MoreArgs=list(tmax.dates=tasmax.dates,tmin.dates=tasmin.dates,prec.dates=pr.dates,
                                  base=c(1951,2010))) ##as.numeric(strsplit(past.int,'-')[[1]])))

      for (k in seq_along(climdex.names)) {
        climdex.values <- lapply(climdex.objects,climdex.names[k])
        ncol <- length(climdex.values[[1]])
        climdex.matrix <- matrix(unlist(climdex.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
        ncvar_put(clim.ncs[[k]],varid=climdex.vars[k],vals=climdex.matrix,
                  start=c(i,1,1),count=c(1,-1,-1))
      }
  }
  lapply(clim.ncs,nc_close)
}

##**************************************************************************************

  climdex.names <- c('climdex.fd',
                     'climdex.su',
                     'climdex.id',
                     'climdex.tr',
                     'climdex.gsl',
                     'climdex.txx',
                     'climdex.tnx',
                     'climdex.txn',
                     'climdex.tnn',
                     'climdex.tn10p',
                     'climdex.tx10p',
                     'climdex.tn90p',
                     'climdex.tx90p',
                     'climdex.wsdi',
                     'climdex.csdi',
                     'climdex.dtr',
                     'climdex.rx1day',
                     'climdex.rx5day',
                     'climdex.sdii',
                     'climdex.r10mm',
                     'climdex.r20mm',
                     'climdex.cdd',
                     'climdex.cwd',
                     'climdex.r95ptot',
                     'climdex.r99ptot',
                     'climdex.prcptot')

##**************************************************************************************
##**************************************************************************************
##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------
##Climdex from the 10km BCCAQ-GCM output
run.anusplin <- function() {

  data.dir <-  '/storage/data/scratch/ssobie/bccaq_gcm_van_whistler_subset/test/'
  write.dir <- '/storage/data/scratch/ssobie/bccaq_gcm_van_whistler_subset/test/rcp85/climdex/'
  
  past.int <- '1951-2000'
  
  climdex.names <- sort(climdex.names)

    gcm <- 'ACCESS1-0'
    rcm <- NULL ##model[2]
    print(gcm)

  first <- lapply(climdex.names,create.climdex.base.files,
                  gcm,rcm=NULL,type=NULL,
                  past.int,
                  data.dir,write.dir)
  
  second <- climdex.for.model(gcm,rcm,interval,type=NULL,
                              climdex.names,
                              past.int,
                              data.dir,write.dir)  
}

run.anusplin()

