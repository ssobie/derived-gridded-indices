##Script to calculate the degree-day indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

library(ncdf4)
library(PCICt)

library(doParallel)
registerDoParallel(cores=3) 

##--------------------------------
##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}

fd <- function(tmin) {
  days <- sum(tmin>0,na.rm=T)
  return((days))
}

s3 <- function(tmax) {
  days <- sum(tmax>30,na.rm=T)
  return((days))
}


gdd<-function(data,fac){tapply(data,fac, dd, tbase=5)}   ##Growing degree days
cdd<-function(data,fac){tapply(data,fac, dd, tbase=18)}  ##Cooling degree days
hdd<-function(data,fac){tapply(-data,fac,dd, tbase=-18)} ##Heating degree days
fdd<-function(data,fac){tapply(-data,fac,dd, tbase=0)} ##Freezing degree days
ffd<-function(data,fac){tapply(data,fac,fd)} ##Frost Free days
s30<-function(data,fac){tapply(data,fac,s3)} ##Days over 30 degC

dd.fxns <- list(gdd=gdd,
                cdd=cdd,   
                hdd=hdd,
                fdd=fdd,
                s30=s30)


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


create.base.files <- function(degree.name,gcm,scenario,type=NULL,
                              past.int,proj.int,new.int,
                              data.dir,write.dir,tmp.dir) {

  ##files <- list.files(path=data.dir,pattern=gcm,full.name=TRUE)
  scen.all <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmax_day')
  ##files.all <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmax_gcm_prism',full.name=TRUE)
  files.all <- scen.all[grep(scenario,scen.all)]
  past.file <- files.all[grep(past.int,files.all)]
  proj.file <- files.all[grep(proj.int,files.all)]

  file.split <- strsplit(past.file,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  ##write.clim.name <- paste(degree.name,'_annual_',gcm,'_',scenario,'_',run,'_',new.int,'.nc',sep='')
  write.clim.name <- paste(degree.name,'_seasonal_',gcm,'_',scenario,'_',run,'_',new.int,'.nc',sep='')
  write.dir <- paste(write.dir,scenario,'/degree_days/',gcm,'/',sep='')

  file.copy(from=paste0(data.dir,gcm,'/',past.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,gcm,'/',proj.file),to=tmp.dir,overwrite=TRUE)

  nc <- nc_open(paste0(tmp.dir,past.file),write=FALSE)
  fnc <- nc_open(paste0(tmp.dir,proj.file),write=FALSE)  
  
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
  full.series <- format(past.origin + (full.values)*86400,'%Y-%m-%d')
  
  years.ix <- grep('*-01-01',full.series)
  years <- full.values[years.ix]
  months.ix <- grep('[0-9]{4}-[0-9]{2}-01',full.series) ##grep('*-*-01',full.series)
  months <- full.values[months.ix]

  seas.ix <- grepl('(-01-01|-04-01|-07-01|-10-01)',full.series)
  seasons <- full.values[seas.ix]
  ##dates <- as.numeric(years)

  dates <- as.numeric(seasons)
  
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

  var.geog <- ncvar_def(degree.name, units='degree days', dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  
  file.nc <- nc_create(paste(tmp.dir,write.clim.name,sep=''), var.geog)
  
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
  ncatt_put(file.nc,varid=degree.name,attname='units',attval='degree days')
  ncatt_put(file.nc,varid=degree.name,attname='_FillValue',attval=-32768)
  ncatt_put(file.nc,varid=degree.name,attname='standard_name',attval=toupper(degree.name))
  ncatt_put(file.nc,varid=degree.name,attname='long_name',attval=toupper(degree.name))

  nc_close(file.nc)
  nc_close(nc)

  if (!file.exists(write.dir))
    dir.create(write.dir,recursive=TRUE)

  file.copy(from=paste0(tmp.dir,write.clim.name),to=write.dir,overwrite=TRUE)
  file.remove(paste0(tmp.dir,past.file))
  file.remove(paste0(tmp.dir,proj.file))

}


##---------------------------------------------------------------

degree.days.for.model <- function(gcm,scenario,interval,type=NULL,
                                  degree.names,
                                  past.int,proj.int,new.int,
                                  data.dir,write.dir,tmp.dir) {
  tasmax.scen <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmax_day')   
  tasmax.files <- tasmax.scen[grep(scenario,tasmax.scen)]
  ###tasmax.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmax_gcm_prism')
  tasmax.past.file <- tasmax.files[grep(past.int,tasmax.files)]
  tasmax.proj.file <- tasmax.files[grep(proj.int,tasmax.files)]
 
  tasmin.scen <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmin_day')
  tasmin.files <- tasmin.scen[grep(scenario,tasmin.scen)]
###  ##tasmin.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='tasmin_gcm_prism')
  tasmin.past.file <- tasmin.files[grep(past.int,tasmin.files)]
  tasmin.proj.file <- tasmin.files[grep(proj.int,tasmin.files)]

  file.split <- strsplit(tasmax.past.file,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  hist.dir <- paste(write.dir,scenario,'/degree_days/',gcm,'/',sep='')    

  ##degree.files <- list.files(path=hist.dir,pattern='annual',full.name=TRUE)
  degree.files <- list.files(path=hist.dir,pattern='seasonal')

  file.copy(from=paste0(data.dir,gcm,'/',tasmax.past.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,gcm,'/',tasmax.proj.file),to=tmp.dir,overwrite=TRUE)

  file.copy(from=paste0(data.dir,gcm,'/',tasmin.past.file),to=tmp.dir,overwrite=TRUE)
  file.copy(from=paste0(data.dir,gcm,'/',tasmin.proj.file),to=tmp.dir,overwrite=TRUE)

  file.copy(from=paste0(hist.dir,degree.files),to=tmp.dir,overwrite=TRUE)

  clim.ncs <- lapply(paste0(tmp.dir,degree.files),nc_open,write=TRUE)

  ##--------------------------------------------------------------

  tasmax.past.nc <- nc_open(paste0(tmp.dir,tasmax.past.file),write=FALSE)
  tasmin.past.nc <- nc_open(paste0(tmp.dir,tasmin.past.file),write=FALSE)

  tasmax.proj.nc <- nc_open(paste0(tmp.dir,tasmax.proj.file),write=FALSE)
  tasmin.proj.nc <- nc_open(paste0(tmp.dir,tasmin.proj.file),write=FALSE)
  
  tn.units <- ncatt_get(tasmax.past.nc,'tasmax')$units
  temp.offset <- 0
  if (tn.units == 'K')
    temp.offset <- 273
  
  ##Attributes to retain
  lon <- ncvar_get(tasmax.past.nc,'lon')
  lat <- ncvar_get(tasmax.past.nc,'lat')  
  n.lon <- length(lon)
  n.lat <- length(lat)

  ##Combine the dates
  time.atts <- ncatt_get(tasmax.past.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                          cal=time.calendar)
  proj.time.atts <- ncatt_get(tasmax.proj.nc,'time')
  proj.time.calendar <- proj.time.atts$calendar
  proj.time.units <- proj.time.atts$units
  
  proj.origin <- as.PCICt(strsplit(proj.time.units, ' ')[[1]][3],
                          cal=proj.time.calendar)
  
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
  yearly.fac <- as.factor(format(tasmax.dates,'%Y'))
  seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")

  seasonal.fac <- get.seasonal.fac(tasmax.dates)
  seas.fac <- seasonal.fac$fac

  ##--------------------------------------------------------------
  ##Compute degree day values and load into newly created netcdf
  for (i in 1:n.lon) {
    print(paste('Lon: ',i,' in ',n.lon,sep=''))
#    for (j in 1:n.lat) {
    tasmax.past.subset <- ncvar_get(tasmax.past.nc,'tasmax',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
    tasmin.past.subset <- ncvar_get(tasmin.past.nc,'tasmin',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
    
    tasmax.proj.subset <- ncvar_get(tasmax.proj.nc,'tasmax',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
    tasmin.proj.subset <- ncvar_get(tasmin.proj.nc,'tasmin',start=c(i,1,1),count=c(1,-1,-1))-temp.offset
        
    tasmax.subset <- cbind(tasmax.past.subset,tasmax.proj.subset)
    tasmin.subset <- cbind(tasmin.past.subset,tasmin.proj.subset)
    
    temp.subset <- (tasmax.subset + tasmin.subset)/2 ##tasmin.subset ##
    ##print('Check for input types')
    
    flag <- is.na(temp.subset[,1])
    print(paste0('Flagged: ',sum(flag),' of ',length(flag)))


    ##Convert temperature subset to list
    temp.list <- lapply(seq_len(nrow(temp.subset)), function(k) temp.subset[k,])
    ##Correct for Dec year in DJF
    avg.fac <- seasonal.fac$fac
    dec.fix <- seasonal.fac$fix
    temp.corr.list <- lapply(temp.list,function(x,y){x[y]},dec.fix)
   

    ##--------------------------------------------------------------
    ##All NA Values
    if (sum(flag) == length(flag)) {
       print('All NA Values')
       temp.list <- vector(mode='list',length=n.lat)
       dlen <- length(levels(seas.fac[[1]])) * length(levels(seas.fac[[2]]))
       degree.values <- lapply(temp.list,function(x){return(as.numeric(rep(NA,dlen)))})
       ncol <- length(degree.values[[1]])
       degree.matrix <- matrix(unlist(degree.values),nrow=n.lat,ncol=ncol,byrow=TRUE)       

       for (k in seq_along(degree.names)) {
          degree.name <- degree.names[k]
          clim.nc <- clim.ncs[[k]]
          ncvar_put(clim.ncs[[k]],varid=degree.names[k],vals=degree.matrix,
                    start=c(i,1,1),count=c(1,-1,-1))
       }

    } else {    ##No or Some NA Values

    ##--------------------------------------------------------------

       if (sum(flag) == 0) {
          print('No NA Values')
          for (k in seq_along(degree.names)) {
             degree.name <- degree.names[k]
             clim.nc <- clim.ncs[[k]]
             fx <- dd.fxns[[degree.names[k]]]

             degree.values <- foreach(
                              temp=temp.corr.list,
                              .export=c('seas.fac','fx')
                              ) %dopar% { 
                                   degree.values <- fx(temp,seas.fac)
                              }
             degree.series <- lapply(degree.values,function(x){as.vector(t(x))})
             ncol <- length(degree.values[[1]])
             degree.matrix <- matrix(unlist(degree.series),nrow=n.lat,ncol=ncol,byrow=TRUE)       
             ncvar_put(clim.ncs[[k]],varid=degree.names[k],vals=degree.matrix,
                       start=c(i,1,1),count=c(1,-1,-1))
          }    

       } else { ##Some NA Values
          print('Some NA Values')
          temp.sub.list <- temp.corr.list[!flag]
          for (k in seq_along(degree.names)) {
             degree.name <- degree.names[k]
             clim.nc <- clim.ncs[[k]]
             fx <- dd.fxns[[degree.names[k]]]

             degree.values <- foreach(
                              temp=temp.sub.list,
                              .export=c('seas.fac','fx')
                              ) %dopar% { 
                                   degree.values <- fx(temp,seas.fac)
                              }
             degree.series <- lapply(degree.values,function(x){as.vector(t(x))})
             all.values <- lapply(temp.list,function(x){return(as.numeric(rep(NA,dlen)))})
             all.values[!flag] <- degree.series
             degree.matrix <- matrix(unlist(all.values),nrow=n.lat,ncol=ncol,byrow=TRUE)       
             ncvar_put(clim.ncs[[k]],varid=degree.names[k],vals=degree.matrix,
                       start=c(i,1,1),count=c(1,-1,-1))
          } ##Degree Loop
       } ##Some NA
    } ##All or none NA
  }##Longitude Loop  
  lapply(clim.ncs,nc_close)



  file.copy(from=paste0(tmp.dir,degree.files),to=hist.dir,overwrite=TRUE)
  file.remove(paste0(tmp.dir,tasmax.past.file))
  file.remove(paste0(tmp.dir,tasmax.proj.file))
  file.remove(paste0(tmp.dir,tasmin.past.file))  
  file.remove(paste0(tmp.dir,tasmin.proj.file))
  file.remove(paste0(tmp.dir,degree.files))

}##Function end

##**************************************************************************************

degree.names <- 's30'
degree.names <- c('cdd',
                  'gdd',
                  'hdd',
                  'fdd')              

  
##  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/'
##  write.dir <- '/home/data/projects/rci/data/stat.downscaling/scaling_comparison/gcm/climdex/'
## 
##
## 
## 

  gcm.list <- c('CanESM2',
                'CCSM4',
                'CNRM-CM5',
                'CSIRO-Mk3-6-0',
                'GFDL-ESM2G',
                'HadGEM2-ES',
                'MIROC5',
                'MPI-ESM-LR',
                'MRI-CGCM3')
##'ACCESS1-0',
  gcm.list <- c('CanESM2',
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
gcm.list <- 'ACCESS1-0'



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
run.bccaq.raw <- function(gcm,scenario) {
##Running BC only versions
  
  data.dir <- paste('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',sep='')
  write.dir <- paste('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',sep='')
  tmp.dir <- '/local_temp/ssobie/dd/'
  if (!file.exists(tmp.dir))
    dir.create(tmp.dir,recursive=T)
 
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'
  
  degree.names <- sort(degree.names)
  print(gcm)

    first <- lapply(degree.names,create.base.files,
                    gcm,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    data.dir,write.dir,tmp.dir)

  second <- degree.days.for.model(gcm,scenario,interval,type='bccaq',
                                  degree.names,
                                  past.int,proj.int,new.int,
                                  data.dir,write.dir,tmp.dir)
}

##-----------------------------------------------------------------------------
##Degree Days from the 800m PRISM adjusted BCCAQ
run.bccaq.prism <- function() {
##Running BC only versions
  
  scenario <- 'rcp85'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'
  region <- 'nanaimo'

  data.dir <- paste('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/',sep='')
  write.dir <- paste('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/',sep='')

   ##Move data to local storage for better I/O
  tmp.dir <- paste('/local_temp/ssobie/',region,'/',sep='')
  if (!file.exists(tmp.dir))
    dir.create(tmp.dir,recursive=TRUE)    
  
  degree.names <- sort(degree.names)
  if ('s30' %in% degree.names) {
    print('Make sure max temp is used for the input')
    ##browser()
  }

  print(gcm.list)
  for (gcm in gcm.list) {
    print(gcm)

    move.to <- paste("rsync -av ",data.dir,gcm,"/*tasmax_gcm_prism* ",tmp.dir,gcm,sep='')
    ##move.to <- paste("rsync -av ",data.dir,gcm,"/tasmax_day* ",tmp.dir,gcm,sep='')
    print(move.to)
    system(move.to)

    print('Create new files')
    first <- lapply(degree.names,create.base.files,
                    gcm,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    tmp.dir,tmp.dir)

    print('Calculate Degree Days')
    second <- degree.days.for.model(gcm,scenario,interval,type='bccaq',
                                    degree.names,
                                    past.int,proj.int,new.int,
                                    tmp.dir,tmp.dir)
    move.back <- paste("rsync -av ",tmp.dir,scenario,"/degree_days/",gcm," ",write.dir,scenario,"/degree_days",sep='')
    print(move.back)
    system(move.back)

    clean.up <- paste("rm ",tmp.dir,gcm,"/*nc" ,sep='')
    print(clean.up)      
    system(clean.up)

    clean.up.dd <- paste("rm ",tmp.dir,scenario,"/degree_days/",gcm,"/*nc" ,sep='')
    print(clean.up.dd)      
    system(clean.up.dd)

  }  
}

##**************************************************************************************

##run.bccaq.prism()

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

##gcm <- 'ACCESS1-0'
##scenario <- 'rcp45'

print(gcm)
print(scenario)
run.bccaq.raw(gcm,scenario)
