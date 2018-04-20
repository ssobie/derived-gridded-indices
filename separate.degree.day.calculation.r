##Script to calculate and write the standard set of derived variables
##for the 800m data

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

degree.days.for.model <- function(degree.names,dd.ncs,lat.ix,n.lon,flag,
                                  tas.list,yearly.fac) {

    ##Degree Day -could be replaced by the function at the beginning
    for (k in seq_along(degree.names)) {
      print(degree.names[k])
      flen <- sum(!flag)
      degree.name <- degree.names[k]
      dd.nc <- dd.ncs[[k]]
      fx <- dd.fxns[[degree.names[k]]]
      sub.list <- tas.list[!flag]

      ftm <- proc.time()      
      sub.values <- foreach(
                         temp=sub.list,
                         .export=c('yearly.fac','fx','dd')
                         ) %do% {
                              degree.values <- fx(temp,yearly.fac)
                         }
       print('Foreach time:')
       print(proc.time()-ftm)

      ncol <- length(sub.values[[1]])
      degree.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(sub.values),nrow=flen,ncol=ncol,byrow=TRUE)
      rm(sub.values)
      degree.matrix[!flag,] <- sub.matrix
      rm(sub.matrix)

      ncvar_put(dd.ncs[[k]],varid=degree.names[k],vals=degree.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
      rm(degree.matrix)       
      gc()
    }##Degree day loop
}


##----------------------------------------------------------------------------------------------

gsl.for.model <- function(ann.name,ann.ncs,lat.ix,n.lon,flag,
                          ann.list,yearly.fac,dates) {
   ##Variables
       flen <- sum(!flag)
       ann.fx <- ann.fxns[[ann.name]]
       sub.list <- ann.list[!flag]
       rm(ann.list)

       ann.avg.values <- foreach(
                         data=sub.list,
                         .export=c('yearly.fac','ann.fx','dates')
                         ) %do% {
                              ann.avg.values <- ann.fx(data,yearly.fac,dates)
                         }
      rm(sub.list)
      ncol <- length(ann.avg.values[[1]])
      ann.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(ann.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
      rm(ann.avg.values)
      ann.avg.matrix[!flag] <- sub.matrix
      rm(sub.matrix)
      ncvar_put(ann.ncs,varid=paste0(ann.name,'ETCCDI'),vals=ann.avg.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
      rm(ann.avg.matrix)
      gc()
}

##----------------------------------------------------------------------------------------------
dtr.for.model <- function(mon.name,mon.ncs,lat.ix,n.lon,flag,
                                   mon.list,monthly.fac) {

   ##Variables
   flen <- sum(!flag)
   mon.fx <- mon.fxns[[mon.name]]
   sub.list <- mon.list[!flag]
   rm(mon.list)

   mon.avg.values <- foreach(
                             data=sub.list,
                             .export=c('monthly.fac','mon.fx')
                             ) %do% {
                                 mon.avg.values <- mon.fx(data,monthly.fac)
                             }
   rm(sub.list)

   ncol <- length(mon.avg.values[[1]])
   mon.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
   sub.matrix <- matrix(unlist(mon.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
   rm(mon.avg.values)
   mon.avg.matrix[!flag,] <- sub.matrix
   rm(sub.matrix)
   ncvar_put(mon.ncs,varid=paste0(mon.name,'ETCCDI'),vals=mon.avg.matrix,
            start=c(1,lat.ix,1),count=c(-1,1,-1))
   rm(mon.avg.matrix)
   gc()
}


##--------------------------------------------------------------
##****************************************************************

if (1==0) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

tmp.dir <- '/local_temp/ssobie/prism/'##tmpdir ## ##tmpdir

gcm <- 'CNRM-CM5'
scenario <- 'rcp85'
run <- 'r1i1p1'
interval <- '1951-2100'
type <- 'gsl'
varname <- 'tas'

##Latitude Bands
###lat.st <- seq(48.1,59.9,0.1)
###lat.en <- seq(48.2,60.0,0.1)

lat.st <- seq(59.6,59.9,0.1)
lat.en <- seq(59.7,60.0,0.1)

len <- length(lat.st)

if (cedar) {
   base.dir <- '/scratch/ssobie/prism/'
} else {
   base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
}

###data.dir <- paste0(base.dir,gcm,'/lat_split/')
data.dir <- paste0('/storage/data/climate/downscale/CMIP5_delivery/lat_split/')


##Move data to local storage for better I/O
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

dir.create(paste0(tmp.dir,scenario),recursive=T)


rtm <- proc.time()
print('Copying')

if (type=='degree_days') {
  template.dir <- paste0(base.dir,gcm,'/template/',scenario,'/degree_days/')
  tmp.dir <- paste0(tmp.dir,scenario,'/')
  dir.create(tmp.dir,recursive=T)  
  ##file.copy(from=template.dir,to=tmp.dir,overwrite=TRUE,recursive=TRUE)
}

if (type=='gsl') {
  template.dir <- paste0(base.dir,gcm,'/template/',scenario,'/climdex')
  template.file <- list.files(path=template.dir,pattern=paste0('gslETCCDI'),full.name=TRUE)
  tmp.dir <- paste0(tmp.dir,scenario,'/climdex/')
  dir.create(tmp.dir,recursive=T)  
  print(template.file)
  ##file.copy(from=template.file,to=tmp.dir,overwrite=TRUE) ##
}

if (type=='dtr') {
  template.dir <- paste0(base.dir,gcm,'/template/',scenario,'/climdex')
  template.file <- list.files(path=template.dir,pattern=paste0('dtrETCCDI'),full.name=TRUE)
  tmp.dir <- paste0(tmp.dir,scenario,'/climdex/')
  dir.create(tmp.dir,recursive=T)  
  print(template.file)
  ##file.copy(from=template.file,to=tmp.dir,overwrite=TRUE) ##
}

print('Move template time')
print(proc.time()-rtm)

##---------------------------------------------------------------------------
  ##Degree Day Files for writing
if (type=='degree_days') {
  print('Degree days opening')
  degree.names <- c('cdd','fdd','gdd','hdd')
  dd.dir <- paste0(tmp.dir,'degree_days/')
  dd.ncs <- vector(mode='list',length=length(degree.names))
  for (d in seq_along(degree.names)) {
    dd.file <- paste0(dd.dir,degree.names[d],'_annual_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
    dd.ncs[[d]] <- nc_open(dd.file,write=TRUE)
  }
  common.lat <- ncvar_get(dd.ncs[[1]],'lat')
  tas <- TRUE
  tasdiff <- FALSE
}

  ##GSL Files for writing
if (type=='gsl') {
  print('GSL opening')
  gsl.dir <- tmp.dir
  gsl.file <- paste0(gsl.dir,'gslETCCDI_ann_BCCAQ2-PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  gsl.ncs <- nc_open(gsl.file,write=TRUE)
  common.lat <- ncvar_get(gsl.ncs,'lat')
  tas <- TRUE
  tasdiff <- FALSE
}

if (type=='dtr') {
  ##DTR Files for writing
  print('DTR opening')
  dtr.dir <- tmp.dir
  dtr.file <- paste0(dtr.dir,'dtrETCCDI_mon_BCCAQ2-PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  dtr.ncs <- nc_open(dtr.file,write=TRUE)
  common.lat <- ncvar_get(dtr.ncs,'lat')
  tas <- FALSE
  tasdiff <- TRUE
}

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

  print('TX opening')
  tasmax.file <- paste0('tasmax_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmax.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
  tasmax.dates <- netcdf.calendar(tasmax.nc)
  yearly.fac <- as.factor(format(tasmax.dates,'%Y'))
  monthly.fac <- as.factor(format(tasmax.dates,'%Y-%m'))

  print('TN Opening')
  tasmin.file <- paste0('tasmin_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmin.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)

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
    if (tas) {
      input.subset <- (tasmax.subset + tasmin.subset)/2
      input.list <- vector(mode='list',length=n.lon)
      input.list <- lapply(seq_len(nrow(input.subset)), function(i) input.subset[i,])
    }
    if (tasdiff) {
       input.subset <- tasmax.subset - tasmin.subset
       input.list <- vector(mode='list',length=n.lon)
       input.list <- lapply(seq_len(nrow(input.subset)), function(i) input.subset[i,])
    }
    flag <- is.na(input.subset[,1])
    rm(input.subset)
    rm(tasmax.subset)
    rm(tasmin.subset)


    ##----------------------------------------------------------
    if (type=='degree_days') {
      ##Degree Day 
      degree.days.for.model(degree.names,dd.ncs,lat.ix,n.lon,flag,
                            input.list,yearly.fac)
    }
    ##----------------------------------------------------------
    ##GSL
    if (type=='gsl') {
      gsl.for.model('gsl',gsl.ncs,lat.ix,n.lon,flag,
                    input.list,yearly.fac,tasmax.dates)
    }
 
    ##----------------------------------------------------------
    ##DTR 
    if (type=='dtr') {
      dtr.for.model('dtr',dtr.ncs,lat.ix,n.lon,flag,
                    input.list,monthly.fac)
    }

    print('Lon loop time')
    print(proc.time()-ltm)

  }##Longitude Loop
  nc_close(tasmax.nc)
  nc_close(tasmin.nc)

  print('Removing lat band files')
  file.remove(paste0(tmp.dir,"/",tasmax.input))
  file.remove(paste0(tmp.dir,"/",tasmin.input))

}##Latitude File Loop

##Move back
write.dir <- paste0(base.dir,gcm)
if (type=='degree_days') {
  for (d in seq_along(degree.names)) {
    nc_close(dd.ncs[[d]]) 
  }
browser()
  file.copy(from=paste0(tmp.dir,scenario,"/degree_days/"),to=paste0(write.dir,'/',scenario,'/'),overwrite=TRUE,recursive=TRUE)

  file.remove(dd.dir)
}


if (type=='gsl') {
  nc_close(gsl.ncs) 
  file.copy(from=gsl.file,to=paste0(write.dir,'/',scenario,'/climdex/'),overwrite=TRUE)
  file.remove(gsl.file)
}
if (type=='dtr') {
  nc_close(dtr.ncs) 
  file.copy(from=dtr.file,to=paste0(write.dir,'/',scenario,'/climdex/'),overwrite=TRUE)
file.remove(dtr.file)
}

print("Total Elapsed Time:")
print(proc.time()-ptm)

