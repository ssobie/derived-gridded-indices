##Script to calculate return periods from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file
##Modified to add the confidence intervals to the resulting netcdf file

library(ncdf4)
library(PCICt)
library(extRemes)
library(ismev)
library(udunits2)
library(zoo)

calc.return.periods <- function(ts.yearly,varname,rperiod) {
  
  if (sum(is.na(ts.yearly)) == length(ts.yearly)) {
    return(NA)
  } else {

    inf.flag <- is.infinite(ts.yearly)
    ts.yearly[inf.flag] <- NA

    na.flag <- is.na(ts.yearly)
    if(sum(na.flag)>0) {
      ts.to.fit <- as.vector(ts.yearly[-which(na.flag)])
    } else {
      ts.to.fit <- as.vector(ts.yearly)
    }
    if (varname=='tasmin') {
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
    rv <- as.numeric(ts.rp)

    if (varname=='tasmin') {
      rv <- -as.numeric(ts.rp) ##$return.level
    }

    return(rv)
  }
}


make.new.netcdf.file <- function(gcm,scenario,varname,rperiod,
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
  
  var.atts <- ncatt_get(nc,varname)
  
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
  varnames <- names(var.atts)
  for (j in 1:length(var.atts))
    ncatt_put(hist.nc,varid=rp.name,attname=varnames[j],attval=var.atts[[j]])

  var.units <- ncatt_get(nc,varname,'units')$value
  ncatt_put(hist.nc,varid=rp.name,attname='units',attval=var.units)   ##'kg m-2 d-1')
  
  nc_close(hist.nc)  
  nc_close(nc)  

}

check.model.outliers <- function(data,varname) {
  rv <- data
  if (varname=='tasmax')
    flags <- which(data > 75)
  if (varname=='tasmin')
    flags <- which(data < -75)
  if (varname=='pr')
    flags <- which(data > 500)
  if (length(flags)!=0)
    rv[flags] <- NA

  return(rv)
}


check.rp.outliers <- function(data,varname) {
  rv <- data
  if (varname=='tasmax')
    flags <- which(data > 75)
  if (varname=='tasmin')
    flags <- which(data < -75)
  if (varname=='pr')
    flags <- which(data > 500)
  if (length(flags)!=0)
    rv[flags] <- 1111

  return(rv)
}


rp.for.model <- function(gcm,scenario,varname,rperiod,
                         var.file,write.file,
                         data.dir,write.dir,
                         interval=NULL) {

  rp.name <- paste('rp.',rperiod,sep='')  

  nc <- nc_open(paste(write.dir,write.file,sep=''),write=TRUE)
  hist.nc <- nc_open(var.file,write=FALSE)
  var.units <- ncatt_get(hist.nc,varname,'units')$value
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

  print(n.lon)
  for (i in 1:n.lon) {
      print(paste0(i,' in ',n.lon))
      var.subset <- ncvar_get(hist.nc,varname,start=c(i,1,st),count=c(1,-1,(en-st+1)))
      if (sum(is.na(var.subset)) != length(data))
        var.subset <- check.model.outliers(var.subset,varname)
      var.list <- list()
      for (j in 1:n.lat) {
        var.list[[j]] <- var.subset[j,]
      }
      if (var.units == 'kg m-2 s-1')
        var.list <- lapply(var.list,ud.convert,var.units,'kg m-2 d-1')

      rp.values <- lapply(var.list,calc.return.periods,varname=varname,rperiod)
      ##print('rp.values')
      
      ##RP
      ##print(max(unlist(rp.values),na.rm=T))
      rp.rp.checked <- lapply(rp.values,check.rp.outliers,varname)
      rp.rp.write <- matrix(unlist(rp.rp.checked),nrow=n.lat,ncol=1,byrow=TRUE)
      ncvar_put(nc,varid=rp.name,vals=rp.rp.write,
                   start=c(i,1,1),count=c(1,n.lat,1))            
  }
  
  nc_close(hist.nc)
  nc_close(nc)
}

##**************************************************************************************
##-------------------------------------------------------------------------


if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

##  gcm <- 'CCSM4'                   
##  varname <- 'tasmin'
##  scenario <- 'rcp85'
##  interval <- '1971-2000'
##  rperiod <- '20'

  tmp.base <- tmpdir ##paste('/local_temp/ssobie/bc/',varname,'/',sep='') ##

  data.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',gcm,'/',scenario,'/annual_extremes/')
  write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',gcm,'/',scenario,'/return_periods/')

  tmp.rp <- paste(tmp.base,scenario,'/return_periods/',sep='')
  tmp.rp <- paste0(tmp.base,paste(gcm,varname,rperiod,interval,sep='_'),'/')

    model <- gcm
    print(model)
    gcm <- model[1]
    rcm <- NULL
    if (!file.exists(tmp.rp))
       dir.create(tmp.rp,recursive=TRUE)
    move.to <- paste("rsync -av ",data.dir,varname,"_annual* ",tmp.rp,sep='')
    print(move.to)
    system(move.to)

    write.rp <- paste(tmp.rp,gcm,'/',sep='')    
    if (!file.exists(write.rp))
      dir.create(write.rp,recursive=TRUE)
    
    var.file <- list.files(path=paste(tmp.rp,'/',sep=''),pattern=paste(varname,'_annual',sep=''),full.name=TRUE)

    ##-------------------------------------------------    
    ##Past File                                                      
    file.split <- strsplit(var.file,'_')[[1]]
    run <- file.split[grep('r*i1p1',file.split)]
    write.hist.name <- paste(varname,'_RP',rperiod,'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc',sep='')

    make.new.netcdf.file(gcm,scenario,varname,rperiod,
                         var.file,write.hist.name,
                         tmp.rp,write.rp)
    print('made new file')
    test <- rp.for.model(gcm,scenario,varname,rperiod,
                         var.file,write.hist.name,
                         tmp.rp,write.rp,interval=interval)
     
  move.back <- paste("rsync -av ",tmp.rp,gcm," ",write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up <- paste("rm ",tmp.rp,"/",varname,"_annual_* " ,sep='')
  print(clean.up)
  system(clean.up)

  clean.up.dd <- paste("rm ",write.rp,varname,"_RP*" ,sep='')
  print(clean.up.dd)
  system(clean.up.dd)



##-------------------------------------------------------------------------


