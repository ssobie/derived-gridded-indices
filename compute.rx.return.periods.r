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
    fx <- max

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


make.new.netcdf.file <- function(gcm,scenario,var.name,rperiod,
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

check.rp.outliers <- function(data,var.name) {
  rv <- data
  if (var.name=='tasmax')
    flags <- which(data > 75)
  if (var.name=='tasmin')
    flags <- which(data < -75)
  if (var.name=='pr')
    flags <- which(data > 500)
  if (var.name=='rx5dayETCCDI'|var.name=='rx2dayETCCDI')
    flags <- which(data > 2000)
  if (length(flags)!=0)
    rv[flags] <- 1111

  return(rv)
}


rp.for.model <- function(gcm,scenario,var.name,rperiod,
                         var.file,write.file,
                         data.dir,write.dir,
                         interval=NULL,months) {

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
  if (length(en)==0) {
    en <- length(var.dates)
  }
  dates.sub <- var.dates[st:en]
  mon.sub <- grep(months,dates.sub)

  yearly.fac <- as.factor(format(dates.sub[mon.sub],'%Y'))

  for (i in 1:n.lon) {
      print(paste0(' Lon i: ',i,' of ',n.lon))
      var.data <- ncvar_get(hist.nc,var.name,start=c(i,1,st),count=c(1,-1,en-st+1))
      var.subset <- var.data[,mon.sub]
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
run.bccaq.prism.rp <- function() {


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

  gcm.list <- 'HadGEM2-ES'
  var.name <- 'rx5dayETCCDI'
  scenario <- 'rcp85'
  past.int <- '1971-2000'
  proj.int <- '2071-2100'
  rperiod <- '10'
  region <- 'van_whistler'
  mon.sub <- '(*-11-02|*-12-02|*-01-02|*-02-02)' ##For RX5Day
  ##mon.sub <- '(*-03-02|*-04-02|*-05-02|*-06-02|*-07-02|*-08-02|*-09-02)' ##For RX2DAY

  data.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/',scenario,'/climdex/')
  write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/',scenario,'/return_periods/')

  tmp.base <- paste('/local_temp/ssobie/',region,'/',sep='')
  tmp.rp <- paste(tmp.base,scenario,'/return_periods/',sep='')

  for (model in gcm.list) {
    print(model)
    gcm <- model[1]
    tmp.dir <- paste(tmp.base,gcm,sep='')
    if (!file.exists(tmp.dir))
       dir.create(tmp.dir,recursive=TRUE)
    move.to <- paste("rsync -av ",data.dir,gcm,"/",var.name,"_mon_* ",tmp.dir,sep='')
    print(move.to)
    system(move.to)

    write.rp <- paste(tmp.rp,gcm,'/',sep='')    
    if (!file.exists(write.rp))
      dir.create(write.rp,recursive=TRUE)
    
    rx.file <- list.files(path=paste(tmp.dir,'/',sep=''),pattern=var.name,full.name=TRUE)

    ##-------------------------------------------------    
    ##Past File                                                      
    file.split <- strsplit(rx.file,'_')[[1]]
    run <- file.split[grep('r*i1p1',file.split)]
    write.hist.name <- paste(var.name,'_RPCI',rperiod,'_BCCAQ_PRISM_',gcm,'_',scenario,'_',run,'_',past.int,'.nc',sep='')

    if (1==0) {
    make.new.netcdf.file(gcm,scenario,var.name,rperiod,
                         rx.file,write.hist.name,
                         tmp.dir,write.rp)
    print('made new file')

    test <- rp.for.model(gcm,scenario,var.name,rperiod,
                         rx.file,write.hist.name,
                         tmp.dir,write.rp,interval=past.int,
                         months=mon.sub)
    } 
    if (1==1) {
    ##-------------------------------------------------
    ##Future File
    write.proj.name <- paste(var.name,'_RPCI',rperiod,'_BCCAQ_PRISM_',gcm,'_',scenario,'_',run,'_',proj.int,'.nc',sep='')
    make.new.netcdf.file(gcm,scenario,var.name,rperiod,
                         rx.file,write.proj.name,
                         tmp.dir,write.rp)
    print('made new file')
    test <- rp.for.model(gcm,scenario,var.name,rperiod,
                         rx.file,write.proj.name,
                         tmp.dir,write.rp,interval=proj.int,
                         months=mon.sub)
  }

  move.back <- paste("rsync -av ",tmp.rp,gcm," ",write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up <- paste("rm ",tmp.dir,"/",var.name,"* " ,sep='')
  print(clean.up)
  system(clean.up)

  clean.up.dd <- paste("rm ",write.rp,var.name,"_RPCI*" ,sep='')
  print(clean.up.dd)
  system(clean.up.dd)


  }##GCM loop
}

##-------------------------------------------------------------------------


run.bccaq.prism.rp() 
