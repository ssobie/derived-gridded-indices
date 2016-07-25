##Script to calculate the degree-day indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables



library(ncdf4)
library(PCICt)
library(climdex.pcic)

##--------------------------------
##Precip as Snow Values
old.pas <- function(pr,tas,fac) {
  Rprof('pas.out')
  tas.pos <- tas
  tas.pos[tas.pos<0]<-0
  frac <- exp(-tas.pos)
  model.data <- pr/10 ##precip in cm
  model.data[tas > 10] <- 0
  warm <- tas <= 10 & tas >=0
  model.data[warm] <- pr[warm]*frac[warm]
  snow <- model.data
  ysnow <- tapply(snow,fac,sum,na.rm=T)
  Rprof(NULL)
  return(round(ysnow,1))
}

pas <- function(pr,tas,fac) {

  tas.pos <- tas
  tas.pos[tas.pos<0]<-0
  model.data <- pr/10 ##precip in cm
  model.data[tas > 0] <- 0
  ysnow <- tapply(model.data,fac,sum,na.rm=T)  
  return(round(ysnow,1))
}

create.base.files <- function(var.name,gcm,scenario,type=NULL,
                              past.int,proj.int,new.int,
                              data.dir,write.dir) {

  ##files <- list.files(path=data.dir,pattern=gcm,full.name=TRUE)
  files.all <- list.files(path=paste(data.dir,scenario,'/',gcm,'/',sep=''),pattern='tas_day',full.name=TRUE)
  ##files.all <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tas_gcm_prism',full.name=TRUE)

  past.file <- files.all[grep(past.int,files.all)]
  proj.file <- files.all[grep(proj.int,files.all)]
  run <- strsplit(past.file,'_')[[1]][9]  

  write.clim.name <- paste(var.name,'_',gcm,'_',scenario,'_',run,'_',new.int,'.nc',sep='')
  write.dir <- paste(write.dir,scenario,'/pas/',gcm,'/',sep='')
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

  proj.time.atts <- ncatt_get(fnc,'time')
  proj.time.calendar <- proj.time.atts$calendar
  proj.time.units <- proj.time.atts$units  
  proj.origin <- as.PCICt(strsplit(proj.time.units, ' ')[[1]][3],
                          cal=proj.time.calendar)
  
  past.values <- ncvar_get(nc,'time')
  proj.values <- ncvar_get(fnc,'time')

  full.values <- seq(past.values[1],tail(proj.values,1),by=1) ##seq(past.values[1],tail(past.values,1),by=1) ##
  full.series <- format(past.origin + (full.values-1)*86400,'%Y-%m-%d')

  years.ix <- grep('*-01-01',full.series)
  years <- full.values[years.ix]
  months.ix <- grep('[0-9]{4}-[0-9]{2}-01',full.series) ##grep('*-*-01',full.series)
  months <- full.values[months.ix]
  
  dates <- as.numeric(years)

    ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  

  lon.atts <- ncatt_get(nc,'lon')
  lat.atts <- ncatt_get(nc,'lat')
  global.atts <- ncatt_get(nc,varid=0)
  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, dates,
                      unlim=TRUE, calendar=time.calendar)

  var.geog <- ncvar_def(var.name, units='cm', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  
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
  ncatt_put(file.nc,varid=var.name,attname='units',attval='cm')
  ncatt_put(file.nc,varid=var.name,attname='_FillValue',attval=-32768)
  ncatt_put(file.nc,varid=var.name,attname='standard_name',attval=toupper(var.name))
  ncatt_put(file.nc,varid=var.name,attname='long_name',attval=toupper(var.name))
  print('Created new file')
  nc_close(file.nc)
  nc_close(nc)
}


##---------------------------------------------------------------

pas.for.model <- function(gcm,scenario,interval,type=NULL,
                          var.names,
                          past.int,proj.int,new.int,
                          data.dir,write.dir) {

  pr.files <- list.files(path=paste(data.dir,scenario,'/',gcm,sep=''),pattern='pr_day',full.name=TRUE)
  pr.past.file <- pr.files[grep(past.int,pr.files)]
  pr.proj.file <- pr.files[grep(proj.int,pr.files)]
  
  tas.files <- list.files(path=paste(data.dir,scenario,'/',gcm,sep=''),pattern='tas_day',full.name=TRUE)
  tas.past.file <- tas.files[grep(past.int,tas.files)]
  tas.proj.file <- tas.files[grep(proj.int,tas.files)]

  run <- strsplit(tas.past.file,'_')[[1]][9]
  hist.dir <- paste(write.dir,scenario,'/pas/',gcm,'/',sep='')    

  pas.files <- list.files(path=hist.dir,pattern='pas',full.name=TRUE)
  clim.ncs <- nc_open(pas.files,write=TRUE)
    
  ##--------------------------------------------------------------
  pr.past.nc <- nc_open(pr.past.file,write=FALSE)
  pr.proj.nc <- nc_open(pr.proj.file,write=FALSE)

  tas.past.nc <- nc_open(tas.past.file,write=FALSE)
  tas.proj.nc <- nc_open(tas.proj.file,write=FALSE)
  
  tx.units <- ncatt_get(tas.past.nc,'tas')$units
  temp.offset <- 0
  if (tx.units == 'K')
    temp.offset <- 273
  
  ##Attributes to retain
  lon <- ncvar_get(tas.past.nc,'lon')
  lat <- ncvar_get(tas.past.nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##Combine the dates
  time.atts <- ncatt_get(tas.past.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                          cal=time.calendar)
  proj.time.atts <- ncatt_get(tas.proj.nc,'time')
  proj.time.calendar <- proj.time.atts$calendar
  proj.time.units <- proj.time.atts$units
  
  proj.origin <- as.PCICt(strsplit(proj.time.units, ' ')[[1]][3],
                          cal=proj.time.calendar)
  
  tas.past.values <- ncvar_get(tas.past.nc,'time')
  tas.proj.values <- ncvar_get(tas.proj.nc,'time')
  tas.dates <- c(past.origin + tas.past.values*86400,
                    proj.origin + tas.proj.values*86400)

  yearly.fac <- as.factor(format(tas.dates,'%Y'))

  ##--------------------------------------------------------------
  ##Compute degree day values and load into newly created netcdf
  for (i in 1:n.lon) {

    print(paste('Lon: ',i,' in ',n.lon,sep=''))

    pr.past.subset <- ncvar_get(pr.past.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))
    pr.proj.subset <- ncvar_get(pr.proj.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))

    tas.past.subset <- ncvar_get(tas.past.nc,'tas',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
    tas.proj.subset <- ncvar_get(tas.proj.nc,'tas',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
    pr.subset <- cbind(pr.past.subset,pr.proj.subset)
    tas.subset <- cbind(tas.past.subset,tas.proj.subset)

    temp.list <- vector(mode='list',length=n.lat)
    pr.list <- vector(mode='list',length=n.lat)
    for (j in 1:n.lat) {
      pr.list[[j]] <- pr.subset[j,]
      temp.list[[j]] <- tas.subset[j,]
    }

    pas.values <- mapply(FUN=var.names,pr.list,temp.list,MoreArgs=list(yearly.fac))
    ncol <- dim(pas.values)[1]
    pas.matrix <- matrix(unlist(pas.values),nrow=n.lat,ncol=ncol,byrow=TRUE)
    ncvar_put(clim.ncs,varid=var.names,vals=pas.matrix,
              start=c(i,1,1),count=c(1,-1,-1))

  }##Longitude Loop  
  nc_close(clim.ncs)

}##Function end

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




##**************************************************************************************
  
##Climdex from the 150km GCM output
run.gcms <- function() {

  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/'
  write.dir <- '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/climdex/'
  
  scenario <- 'rcp45' ##'sresa2' ##'rcp45'
  past.int <- '1971-2000'
  proj.int <- '2041-2070'
  new.int <- '1971-2070'
  
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
##Climdex from the 10km BCCAQ-GCM output
run.bccaq.raw <- function() {
##Running BC only versions
  ptm <- proc.time()  
  ##data.dir <-  '/home/data/climate/downscale/CMIP5/BCCAQ/'
  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_bc_subset/'
  write.dir <- '/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm/'
  
  scenario <- 'rcp85'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'
  
  var.names <- 'pas'
  ##gcm.list <- 'MRI-CGCM3'
  for (gcm in gcm.list) {
    print(gcm)

    first <- lapply(var.names,create.base.files,
                    gcm,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)

    second <- pas.for.model(gcm,scenario,interval,type='bccaq',
                                    var.names,
                                    past.int,proj.int,new.int,
                                    data.dir,write.dir)    
  }
  print('Elapsed time')
  print(proc.time() - ptm)
        
}

##-----------------------------------------------------------------------------
##Degree Days from the 800m PRISM adjusted BCCAQ
run.bccaq.prism <- function() {
##Running BC only versions
  
  ##data.dir <-  '/home/data/climate/downscale/CMIP5/BCCAQ/'

  data.dir <- '/home/data/scratch/ssobie/bccaq_gcm_okanagan_subset/'
  write.dir <- '/home/data/scratch/ssobie/bccaq_gcm_okanagan_subset/'
  
  scenario <- 'rcp85'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'

  gcm.list <- 'MRI-CGCM3'
  var.names <- 'pas'
  var.names <- sort(var.names)

  for (gcm in gcm.list) {
    print(gcm)

    first <- lapply(var.names,create.base.files,
                    gcm,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)

    second <- pas.for.model(gcm,scenario,interval,type='bccaq',
                                    var.names,
                                    past.int,proj.int,new.int,
                                    data.dir,write.dir)
  }  
}

##**************************************************************************************

run.bccaq.raw()
