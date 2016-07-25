##Script to calculate the climdex indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

library(ncdf4)
library(PCICt)
library(climdex.pcic)

climdex.r95days <- function(clim.object) {

  years <- as.Date(paste(unique(format(clim.object@dates,'%Y')),'-01-16',sep=''))
  na.rv <- rep(NA,length(years))
  precip <- clim.object@data$prec
  if (sum(is.na(precip))==length(precip))
    return(na.rv)
  q.obj <- get.outofbase.quantiles(prec=clim.object@data$prec,prec.dates=clim.object@dates,
                                                 base=c(1971,2000))
  q.95 <- number.days.op.threshold(clim.object@data$prec,date.factor=clim.object@date.factors$annual, q.obj$prec[1], op = ">")
  return(q.95)
}

climdex.r99days <- function(clim.object) {

  years <- as.Date(paste(unique(format(clim.object@dates,'%Y')),'-01-16',sep=''))
  na.rv <- rep(NA,length(years))
  precip <- clim.object@data$prec
  if (sum(is.na(precip))==length(precip))
    return(na.rv)
  q.obj <- get.outofbase.quantiles(prec=clim.object@data$prec,prec.dates=clim.object@dates,
                                                 base=c(1971,2000))
  q.99 <- number.days.op.threshold(clim.object@data$prec,date.factor=clim.object@date.factors$annual, q.obj$prec[1], op = ">")
  return(q.99)
}



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
                        climdex.r95days=c('r95daysETCCDI','Ann','days'),
                        climdex.r99days=c('r99daysETCCDI','Ann','days'))
  rv <- climdex.names[[climdex.name]]
  return(rv)
}


get.climdex.value <- function(climdex.var,climdex.object) {
  
  rv <- eval(parse(text=paste(climdex.var,'(climdex.object)',sep='')))
  return(rv)
}


create.climdex.base.files <- function(climdex.name,gcm,rcm=NULL,scenario,type=NULL,
                                      past.int,proj.int,new.int,
                                      data.dir,write.dir) {
  print(climdex.name)
  climdex.info <- get.climdex.info(climdex.name)
  climdex.var <- climdex.info[1]
  climdex.calendar <- climdex.info[2]
  climdex.units <- climdex.info[3]

  if (is.null(rcm)) {
    if (type=='prism') {
      files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='pr_gcm_prism',full.name=TRUE)
      past.file <- files[grep(past.int,files)]
      run <- strsplit(past.file,'_')[[1]][10]
      proj.file <- files[grep(proj.int,files)]     
      
      write.clim.name <- paste(climdex.var,'_',tolower(climdex.calendar),'_BCCAQ-PRISM_',gcm,'_',scenario,'_',run,'_',new.int,'.nc',sep='')
      write.dir <- paste(write.dir,gcm,'/',sep='')
      
    } else {
      files <- list.files(path=paste(data.dir,gcm,'/',sep=''),pattern='pr_day',full.name=TRUE)
##      files <- list.files(path=paste(data.dir,scenario,'/',gcm,'/',sep=''),pattern='pr_day',full.name=TRUE)
      past.file <- files[grep(past.int,files)]

      run <- strsplit(past.file,'_')[[1]][10]
      proj.file <- files[grep(proj.int,files)]     
      
      write.clim.name <- paste(climdex.var,'_',tolower(climdex.calendar),'_',gcm,'_',scenario,'_',run,'_',new.int,'.nc',sep='')
      write.dir <- paste(write.dir,gcm,'/',sep='')
    }
    
    if (!file.exists(write.dir))
      dir.create(write.dir,recursive=TRUE)
    
  } else {
    if (type=='bccaq') {
      read.dir <- paste(data.dir,gcm,'_',rcm,'/',sep='')
      past.file <- paste(read.dir,'pr_day_BCCAQ_',gcm,'_',rcm,'_1971-2000.nc',sep='')
      proj.file <- paste(read.dir,'pr_day_BCCAQ_',gcm,'_',rcm,'_2041-2070.nc',sep='')     
    } 
    if (type=='rcm') {
      read.dir <- paste(data.dir,tolower(gcm),'.',tolower(rcm),'/',sep='')
      past.file <- paste(read.dir,'pr_day_',tolower(gcm),'_',tolower(rcm),'_20c3m_rcm_scale_1971-2000.nc',sep='')
      proj.file <- paste(read.dir,'pr_day_',tolower(gcm),'_',tolower(rcm),'_sresa2_rcm_scale_2041-2070.nc',sep='')     
    }
    write.clim.name <- paste(climdex.var,'_',tolower(climdex.calendar),'_',gcm,'_',rcm,'_',new.int,'.nc',sep='')
    write.dir <- paste(write.dir,gcm,'_',rcm,'/',sep='')
  }

  if (!file.exists(write.dir))
    dir.create(write.dir,recursive=TRUE)

  nc <- nc_open(past.file,write=FALSE)
  fnc <- nc_open(proj.file,write=FALSE)
  
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  past.values <- ncvar_get(nc,'time')

  proj.time.atts <- ncatt_get(fnc,'time')
  proj.time.calendar <- proj.time.atts$calendar
  proj.time.units <- proj.time.atts$units
  
  proj.origin <- as.PCICt(strsplit(proj.time.units, ' ')[[1]][3],
                          cal=proj.time.calendar)
  proj.values <- ncvar_get(fnc,'time')
  
  ##Original time method
  full.values <- seq(past.values[1],tail(proj.values,1),by=1)
  full.series <- format(past.origin + (full.values-1)*86400,'%Y-%m-%d')
    
##  proj.origin <- as.PCICt('2038-01-01',cal=time.calendar)



  ##For GCM Scale from RCM (RCM aggregated to 150km)  
  ##full.values <- seq(past.values[1],tail(proj.values,1),by=1)
  ##proj.series <- format(proj.origin + proj.values*86400,'%Y-%m-%d')

  ##if (as.PCICt('2070-12-30',cal=proj.time.calendar) > as.PCICt(tail(proj.series,1),cal=proj.time.calendar)) {
  ##  fill <- as.numeric(as.PCICt('2070-12-30',cal=proj.time.calendar) - as.PCICt(tail(proj.series,1),cal=proj.time.calendar))/86400
  ##  full.values <- c(full.values,seq(from=tail(proj.values,1)+1,by=1,length.out=fill))
  ##  full.values <- (full.values-1)
  ##  full.series <- format(past.origin + full.values*86400,'%Y-%m-%d')    
  ##}

  ##For RCM or GCM scale from RCM driven output (RCMs to a common grid)
  if (1==0) {
  ##if (type=='rcm') {
    past.series <- format(past.origin + past.values*86400,'%Y-%m-%d')
    proj.series <- format(proj.origin + proj.values*86400,'%Y-%m-%d')
    ##proj.series <- format(proj.origin + floor(past.values)*86400,'%Y-%m-%d')

    time.diff <- as.numeric((proj.origin - past.origin)/86400)
    full.values <- seq(from=head(past.values,1),by=1,length.out=time.diff+length(proj.values))
    
    if (climdex.calendar=='Mon') {
      fill <- as.numeric(as.PCICt('2070-12-30',cal=proj.time.calendar) - as.PCICt(tail(proj.series,1),cal=proj.time.calendar))/86400
      full.values <- seq(from=head(past.values,1),by=1,length.out=time.diff+length(proj.values)+fill)
    }
    full.series <- format(past.origin + full.values*86400,'%Y-%m-%d')
  }
  
  years.ix <- grep('*-01-01',full.series)
  years <- full.values[years.ix]
  months.ix <- grep('[0-9]{4}-[0-9]{2}-01',full.series) ##grep('*-*-01',full.series)
  months <- full.values[months.ix]

  dates <- switch(climdex.calendar,
                  Ann=years,
                  Mon=months)
  dates <- as.numeric(dates)

  
  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  if (type=='rcm')
    lon <- lon - 360
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

climdex.for.model <- function(gcm,rcm=NULL,scenario,interval,type=NULL,
                              climdex.names,
                              past.int,proj.int,new.int,
                              data.dir,write.dir) {

  if (is.null(rcm)) {
    ##pr.files <- list.files(path=paste(data.dir,scenario,'/',gcm,'/',sep=''),pattern='pr_day',full.name=TRUE)
    pr.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='pr_gcm_prism',full.name=TRUE)
    pr.past.file <- pr.files[grep(past.int,pr.files)]
    pr.proj.file <- pr.files[grep(proj.int,pr.files)] 

    run <- strsplit(pr.past.file,'_')[[1]][10] ##[9]

    ##tasmax.files <- list.files(path=paste(data.dir,scenario,'/',gcm,'/',sep=''),pattern='tasmax_day',full.name=TRUE)
    tasmax.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmax_gcm_prism',full.name=TRUE)
    tasmax.past.file <- tasmax.files[grep(past.int,tasmax.files)]
    tasmax.proj.file <- tasmax.files[grep(proj.int,tasmax.files)]  

    ##tasmin.files <- list.files(path=paste(data.dir,scenario,'/',gcm,'/',sep=''),pattern='tasmin_day',full.name=TRUE)
    tasmin.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmin_gcm_prism',full.name=TRUE)
    tasmin.past.file <- tasmin.files[grep(past.int,tasmin.files)]
    tasmin.proj.file <- tasmin.files[grep(proj.int,tasmin.files)]  
       
    hist.dir <- paste(write.dir,gcm,'/',sep='')

    ##clim.all.files <- list.files(path=hist.dir,pattern='r9.daysETCCDI',full.name=TRUE)

    clim.all.files <- list.files(path=hist.dir,pattern='ETCCDI',full.name=TRUE)
    clim.files <- clim.all.files[grep(scenario,clim.all.files)]
    clim.ncs <- lapply(clim.files,nc_open,write=TRUE)

   } else {
     if (type=='bccaq') {
       read.dir <- paste(data.dir,gcm,'_',rcm,'/',sep='')
       pr.past.file <- paste(read.dir,'pr_day_BCCAQ_',gcm,'_',rcm,'_1971-2000.nc',sep='')
       pr.proj.file <- paste(read.dir,'pr_day_BCCAQ_',gcm,'_',rcm,'_2041-2070.nc',sep='')     
       tasmax.past.file <- paste(read.dir,'tasmax_day_BCCAQ_',gcm,'_',rcm,'_1971-2000.nc',sep='')
       tasmax.proj.file <- paste(read.dir,'tasmax_day_BCCAQ_',gcm,'_',rcm,'_2041-2070.nc',sep='')
       tasmin.past.file <- paste(read.dir,'tasmin_day_BCCAQ_',gcm,'_',rcm,'_1971-2000.nc',sep='')     
       tasmin.proj.file <- paste(read.dir,'tasmin_day_BCCAQ_',gcm,'_',rcm,'_2041-2070.nc',sep='')
     }
     if (type=='rcm') {
       read.dir <- paste(data.dir,tolower(gcm),'.',tolower(rcm),'/',sep='')
       pr.past.file <- paste(read.dir,'pr_day_',tolower(gcm),'_',tolower(rcm),'_20c3m_rcm_scale_1971-2000.nc',sep='')
       pr.proj.file <- paste(read.dir,'pr_day_',tolower(gcm),'_',tolower(rcm),'_sresa2_rcm_scale_2041-2070.nc',sep='')     
       tasmax.past.file <- paste(read.dir,'tasmax_day_',tolower(gcm),'_',tolower(rcm),'_20c3m_rcm_scale_1971-2000.nc',sep='')
       tasmax.proj.file <- paste(read.dir,'tasmax_day_',tolower(gcm),'_',tolower(rcm),'_sresa2_rcm_scale_2041-2070.nc',sep='')
       tasmin.past.file <- paste(read.dir,'tasmin_day_',tolower(gcm),'_',tolower(rcm),'_20c3m_rcm_scale_1971-2000.nc',sep='')     
       tasmin.proj.file <- paste(read.dir,'tasmin_day_',tolower(gcm),'_',tolower(rcm),'_sresa2_rcm_scale_2041-2070.nc',sep='')
     }

     hist.dir <- paste(write.dir,gcm,'_',rcm,'/',sep='')
     clim.files <- list.files(path=hist.dir,pattern='ETCCDI',full.name=TRUE) ##list.files(path=hist.dir,pattern='r95days',full.name=TRUE)##
     print(clim.files)
     clim.ncs <- lapply(clim.files,nc_open,write=TRUE)     
  }

  ##--------------------------------------------------------------
  climdex.vars <- sort(unlist(lapply(lapply(climdex.names,get.climdex.info),function(x){return(x[1])})))

  pr.past.nc <- nc_open(pr.past.file,write=FALSE)
  tasmax.past.nc <- nc_open(tasmax.past.file,write=FALSE)
  tasmin.past.nc <- nc_open(tasmin.past.file,write=FALSE)

  pr.proj.nc <- nc_open(pr.proj.file,write=FALSE)
  tasmax.proj.nc <- nc_open(tasmax.proj.file,write=FALSE)
  tasmin.proj.nc <- nc_open(tasmin.proj.file,write=FALSE)

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
  proj.time.atts <- ncatt_get(pr.proj.nc,'time')
  proj.time.calendar <- proj.time.atts$calendar
  proj.time.units <- proj.time.atts$units
  
  proj.origin <- as.PCICt(strsplit(proj.time.units, ' ')[[1]][3],
                           cal=proj.time.calendar)
  
  pr.past.values <- ncvar_get(pr.past.nc,'time')
  pr.proj.values <- ncvar_get(pr.proj.nc,'time')
  pr.time.values <- c(pr.past.values,pr.proj.values)
  pr.dates <- c(past.origin + pr.past.values*86400,
                proj.origin + pr.proj.values*86400)

  tasmax.past.values <- ncvar_get(tasmax.past.nc,'time')
  tasmax.proj.values <- ncvar_get(tasmax.proj.nc,'time')
  tasmax.time.values <- c(tasmax.past.values,tasmax.proj.values)
  tasmax.dates <- c(past.origin + tasmax.past.values*86400,
                    proj.origin + tasmax.proj.values*86400)

  tasmin.past.values <- ncvar_get(tasmin.past.nc,'time')
  tasmin.proj.values <- ncvar_get(tasmin.proj.nc,'time')
  tasmin.time.values <- c(tasmin.past.values,tasmin.proj.values)
  tasmin.dates <- c(past.origin + tasmin.past.values*86400,
                    proj.origin + tasmin.proj.values*86400)

  ##--------------------------------------------------------------
  ##Compute climdex values and load into newly created climdex netcdf
  for (i in 1:n.lon) {
    print(paste('Lon: ',i,' in ',n.lon,sep=''))
#    for (j in 1:n.lat) {
      pr.past.subset <- ncvar_get(pr.past.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))*pr.scale
      tasmax.past.subset <- ncvar_get(tasmax.past.nc,'tasmax',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
      tasmin.past.subset <- ncvar_get(tasmin.past.nc,'tasmin',start=c(i,1,1),count=c(1,-1,-1))-temp.offset

      pr.proj.subset <- ncvar_get(pr.proj.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))*pr.scale
      tasmax.proj.subset <- ncvar_get(tasmax.proj.nc,'tasmax',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
      tasmin.proj.subset <- ncvar_get(tasmin.proj.nc,'tasmin',start=c(i,1,1),count=c(1,-1,-1))-temp.offset

      pr.subset <- cbind(pr.past.subset,pr.proj.subset)
      tasmax.subset <- cbind(tasmax.past.subset,tasmax.proj.subset)
      tasmin.subset <- cbind(tasmin.past.subset,tasmin.proj.subset)

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
                                  base=as.numeric(strsplit(past.int,'-')[[1]])))
    
##      climdex.object <- climdexInput.raw(tasmax.subset,tasmin.subset,pr.subset,
##                                         tasmax.dates,tasmin.dates,pr.dates,
##                                         base=as.numeric(strsplit(past.int,'-')[[1]]))

      for (k in seq_along(climdex.names)) {
        climdex.values <- lapply(climdex.objects,climdex.names[k])
        ncol <- length(climdex.values[[1]])
        climdex.matrix <- matrix(unlist(climdex.values),nrow=n.lat,ncol=ncol,byrow=TRUE)        
        ncvar_put(clim.ncs[[k]],varid=climdex.vars[k],vals=climdex.matrix,
                  start=c(i,1,1),count=c(1,-1,-1))
        print(climdex.names[k])
      }

##      climdex.values <- lapply(climdex.names,get.climdex.value,climdex.object)
##      names(climdex.values) <- climdex.names

##      mapply(ncvar_put,nc=clim.ncs,varid=climdex.vars,vals=climdex.values,
##             MoreArgs=list(start=c(i,j,1),count=c(1,1,-1)))

#    }

  }
  
  lapply(clim.ncs,nc_close)
  
  ##--------------------------------------------------------------

}

##**************************************************************************************

##climdex.names <- c('climdex.tn10p',
##                   'climdex.tx10p',
##                   'climdex.tn90p',
##                   'climdex.tx90p',
##                   'climdex.r95ptot',
##                   'climdex.r99ptot',
##                   'climdex.prcptot')

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
                     'climdex.prcptot',
                     'climdex.r95days',
                     'climdex.r99days')

#  climdex.names <- c('climdex.fd',
#                     'climdex.txx',
#                     'climdex.tnn',
#                     'climdex.tx90p',
#                     'climdex.tn10p',
#                     'climdex.r95ptot',
#                     'climdex.r99ptot',
#                     'climdex.prcptot')

#  climdex.names <- c('climdex.fd')

                   
                   ##'climdex.sdii',
                   ##'climdex.r10mm',
                   ##'climdex.r20mm',
                   ##'climdex.cdd',
                   ##'climdex.cwd')


##    climdex.names <- c('climdex.r95days','climdex.r99days')

##  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/'
##  write.dir <- '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/climdex/'
## 'ACCESS1-0',

## 
## 
## 

##gcm.list <- c('ACCESS1-0',
##              'CanESM2',
##              'CCSM4',
##              'CNRM-CM5',
##              'CSIRO-Mk3-6-0',
##              'GFDL-ESM2G')

gcm.list <- c('HadGEM2-CC',              
              'HadGEM2-ES',
              'inmcm4',
              'MIROC5',
              'MPI-ESM-LR', 
              'MRI-CGCM3')


## 
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

##**************************************************************************************
  
##Climdex from the 150km GCM output
run.gcms <- function() {
  
  scenario <- 'rcp85' ##'sresa2' ##'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  new.int <- '1971-2070'

  data.dir <-  paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/',sep='')
  write.dir <- paste('/home/data/scratch/ssobie/gcm/climdex/',scenario,'/',sep='')
    ##paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/climdex/latest/',sep='')
  
  climdex.names <- sort(climdex.names)
  for (model in gcm.list) {
    gcm <- model[1]
    rcm <- NULL ##model[2]
    
    first <- lapply(climdex.names,create.climdex.base.files,
                    gcm,rcm=NULL,scenario,type='gcm',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)
    
    second <- climdex.for.model(gcm,rcm,scenario,interval,type='gcm',
                                climdex.names,
                                past.int,proj.int,new.int,
                                data.dir,write.dir) 
  }  
}

##Climdex from the 150km GCM output
run.miroc4h <- function() {
  
  scenario <- 'rcp45' ##'sresa2' ##'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2006-2035'
  new.int <- '1971-2035'

  data.dir <-  paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm_bc_subset/',sep='')
  write.dir <- paste('/home/data/scratch/ssobie/gcm/climdex/',scenario,'/',sep='')
    ##paste('/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/',scenario,'/climdex/latest/',sep='')
  
  climdex.names <- sort(climdex.names)
  gcm <- 'MIROC4h'
  rcm <- NULL ##model[2]
    
  first <- lapply(climdex.names,create.climdex.base.files,
                  gcm,rcm=NULL,scenario,type='gcm',
                  past.int,proj.int,new.int,
                  data.dir,write.dir)
  
  second <- climdex.for.model(gcm,rcm,scenario,interval,type='gcm',
                              climdex.names,
                              past.int,proj.int,new.int,
                              data.dir,write.dir) 
}

##**************************************************************************************
##-----------------------------------------------------------------------------
##Climdex from the 150km GCM output
run.bccaq.gcm.scale <- function() {

  data.dir <-  '/home/data/scratch/ssobie/bccaq_gcm/gcm/'
  write.dir <- '/home/data/scratch/ssobie/bccaq_gcm/climdex/'
  
  scenario <- 'rcp45' ##'sresa2' ##'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  new.int <- '1971-2070'
  
  climdex.names <- sort(climdex.names)

  for (model in gcm.list) {
    gcm <- model[1]
    rcm <- NULL ##model[2]
    print(gcm)
    
    first <- lapply(climdex.names,create.climdex.base.files,
                    gcm,rcm=NULL,scenario,type='gcm',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)
    
    second <- climdex.for.model(gcm,rcm,scenario,interval,type='gcm',
                                climdex.names,
                                past.int,proj.int,new.int,
                                data.dir,write.dir) 
  }  
}

##-----------------------------------------------------------------------------
##Climdex from the 150km RCM output
run.bccaq.rcm.scale <- function() {

  data.dir <-  '/home/data/scratch/ssobie/bccaq_rcm/rcm/' ##'/home/data/projects/rci/data/stat.downscaling/scaling_comparison/bccaq_gcm/gcm/'
  write.dir <- '/home/data/scratch/ssobie/bccaq_rcm/rcm/climdex_from_rcm_scale/' ##'/home/data/projects/rci/data/stat.downscaling/scaling_comparison/bccaq_gcm/gcm/climdex/'
  
  scenario <- 'sresa2'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  new.int <- '1971-2070'
  
  climdex.names <- sort(climdex.names)

  for (model in rcm.list) {
    gcm <- model[1]
    rcm <- model[2]
    print(gcm)
    print(rcm)
    first <- lapply(climdex.names,create.climdex.base.files,
                    gcm,rcm,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)
    
    second <- climdex.for.model(gcm,rcm,scenario,interval,type='bccaq',
                                climdex.names,
                                past.int,proj.int,new.int,
                                data.dir,write.dir) 
  }  
}

##-----------------------------------------------------------------------------
run.rcm.rcm.scale <- function() {

  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/rcm/'
  write.dir <- '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/rcm/climdex/'
  
  scenario <- 'sresa2' ##'sresa2' ##'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  new.int <- '1971-2070'
  
  climdex.names <- sort(climdex.names)
  rcm.list <- list(c('CCSM','WRFG'))
  
  for (model in rcm.list) {
    gcm <- model[1]
    rcm <- model[2]
    print(gcm)
    print(rcm)
    
    first <- lapply(climdex.names,create.climdex.base.files,
                    gcm,rcm,scenario,type='rcm',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)
    
    second <- climdex.for.model(gcm,rcm,scenario,interval,type='rcm',
                                climdex.names,
                                past.int,proj.int,new.int,
                                data.dir,write.dir) 
  }  
}

##-----------------------------------------------------------------------------
run.rcm.gcm.scale <- function() {

  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/rcm/gcm/'  ##'/home/data/scratch/ssobie/rcm/gcm/'
  write.dir <- '/home/data/scratch/ssobie/rcm/gcm/climdex/'
  
  scenario <- 'sresa2' ##'sresa2' ##'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  new.int <- '1971-2070'
  
  climdex.names <- sort(climdex.names)

  for (model in rcm.list) {
    gcm <- model[1]
    rcm <- model[2]
    print(gcm)
    print(rcm)
    
    first <- lapply(climdex.names,create.climdex.base.files,
                    gcm,rcm,scenario,type='rcm',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)

    second <- climdex.for.model(gcm,rcm,scenario,interval,type='rcm',
                                climdex.names,
                                past.int,proj.int,new.int,
                                data.dir,write.dir) 

  }  
}

##-----------------------------------------------------------------------------
##Climdex from the 10km BCCAQ-RCM output
run.bccaq.raw <- function() {
  
  scenario <- 'rcp85'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'
  
  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_bc_subset/' ##Scenario is pasted in above
  ##write.dir <- paste('/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm/',scenario,'/climdex/',sep='')
   write.dir <- paste('/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm/rcp85/climdex/full/',sep='')
##  write.dir <- '/home/data/climate/downscale/CMIP5/BCCAQ/climdex/bc_subset/' ##For rcp85  
  ##rcm.list <- list(c('CCSM','WRFG'))
  gcm.list <- 'CanESM2'
  climdex.names <- sort(climdex.names)
  for (model in gcm.list) {
    gcm <- model##[1]
    ##rcm <- model[2]
    print(gcm)
    ##print(paste(gcm,rcm,sep='-'))

    first <- lapply(climdex.names,create.climdex.base.files,
                    gcm,rcm=NULL,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)
    
    second <- climdex.for.model(gcm,rcm=NULL,scenario,interval,type='bccaq',
                                climdex.names,
                                past.int,proj.int,new.int,
                                data.dir,write.dir)
##    browser()
  }  
}

##-----------------------------------------------------------------------------
##Climdex from the 800 BCCAQ-PRISM output
run.bccaq.prism <- function() {
  
  scenario <- 'rcp85'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'
  
  data.dir <-  '/storage/data/scratch/ssobie/bccaq_gcm_van_whistler_subset/' ##Scenario is pasted in above
  write.dir <- paste('/storage/data/scratch/ssobie/bccaq_gcm_van_whistler_subset/rcp85/climdex/',sep='')

  climdex.names <- sort(climdex.names)
  for (model in gcm.list) {
    gcm <- model
    print(gcm)

    first <- lapply(climdex.names,create.climdex.base.files,
                    gcm,rcm=NULL,scenario,type='prism',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)

    second <- climdex.for.model(gcm,rcm=NULL,scenario,interval,type='prism',
                                climdex.names,
                                past.int,proj.int,new.int,
                                data.dir,write.dir)
    
  }  
}



##**************************************************************************************
##Climdex from the 150km GCM output

  ##  data.dir <- '/home/data/climate/downscale/CMIP5/BCCAQ/'
  ##  write.dir <- '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/bccaq_gcm_bc_subset/'

##Rprof('climdex.bccaq.raw.profile.out')
##run.bccaq.rcm.scale()
##run.rcm.rcm.scale()
##Rprof(NULL)
run.bccaq.prism()
