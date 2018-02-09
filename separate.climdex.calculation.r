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
library(doParallel)
registerDoParallel(cores=2)

##----------------------------------------------------------------------------------------------
climdex.extremes.for.model <- function(clim.name,clim.nc,lat.ix,n.lon,flag,
                                      tasmax.list,tasmax.dates,
                                      tasmin.list,tasmin.dates,
                                      pr.list,pr.dates,base) {

   flen <- sum(!flag)
   tasmax.sub <- tasmax.list[!flag]     
   rm(tasmax.list)
   tasmin.sub <- tasmin.list[!flag]
   rm(tasmin.list)
   pr.sub <- pr.list[!flag]
   rm(pr.list)
   climdex.objects <- foreach(
                             tasmax=tasmax.sub,
                             tasmin=tasmin.sub,
                             pr=pr.sub,
                             .export=c('climdexInput.raw','tasmax.dates','tasmin.dates','pr.dates','base')
                             ) %dopar% {
                                objects <- climdexInput.raw(tmax=tasmax,tmin=tasmin,prec=pr,
                                           tmax.dates=tasmax.dates,tmin.dates=tasmin.dates,prec.dates=pr.dates,base=base)
                             }
   
   rm(tasmax.sub)
   rm(tasmin.sub)
   rm(pr.sub)
   
   clim.fx <- clim.fxns[[clim.name]]
   climdex.info <- get.climdex.info(clim.name)

   climdex.values <- foreach(obj=climdex.objects,
                             .export='clim.fx'
                             ) %dopar% {
                             climdex.values <- clim.fx(obj)
                             }
   rm(climdex.objects)
   ncol <- length(climdex.values[[1]])
   climdex.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
   sub.matrix <- matrix(unlist(climdex.values),nrow=flen,ncol=ncol,byrow=TRUE)
   rm(climdex.values)
   climdex.matrix[!flag,] <- sub.matrix
   rm(sub.matrix)
   ncvar_put(clim.nc,varid=climdex.info[1],vals=climdex.matrix,
             start=c(1,lat.ix,1),count=c(-1,1,-1))
   rm(climdex.matrix)                  
   gc()
}


##--------------------------------------------------------------
##****************************************************************

if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

tmp.dir <- tmpdir ##'/local_temp/ssobie/prism/' ##tmpdir

##gcm <- 'CNRM-CM5'
##scenario <- 'rcp85'
##run <- 'r1i1p1'
##interval <- '1951-2100'
##type <- 'climdex'
##varname <- 'fd'

##Latitude Bands
lat.st <- seq(48.1,59.9,0.1) ##format(seq(48.0,59.9,0.1),nsmall=1)
lat.en <- seq(48.2,60.0,0.1) ##format(seq(48.1,60.0,0.1),nsmall=1)

len <- length(lat.st)

if (cedar) {
   base.dir <- '/scratch/ssobie/prism/'
} else {
   base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
}

data.dir <- paste0(base.dir,gcm,'/lat_split/')
template.dir <- paste0(base.dir,gcm,'/template/',scenario,'/',type,'/')
template.file <- list.files(path=template.dir,pattern=paste0(varname,'ETCCDI'),full.name=TRUE)

##Move data to local storage for better I/O
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}
dir.create(paste0(tmp.dir,scenario,'/',type),recursive=T)

rtm <- proc.time()
print('Copying')
print(template.file)
##file.copy(from=template.dir,to=paste0(tmp.dir,scenario,'/'),overwrite=TRUE,recursive=TRUE) ##
file.copy(from=template.file,to=paste0(tmp.dir,scenario,'/',type,'/'),overwrite=TRUE,recursive=TRUE) ##

print('Move template time')
print(proc.time()-rtm)

##---------------------------------------------------------------------------
##Climdex Files for writing
  print('Climdex opening')

clim.name <- paste0('climdex.',varname)

  climdex.dir <- paste0(tmp.dir,'/',scenario,'/climdex/')
  climdex.info <- get.climdex.info(clim.name)
  climdex.var <- climdex.info[1]
  climdex.calendar <- climdex.info[2]
  clim.file <- paste0(climdex.dir,climdex.var,'_',tolower(climdex.calendar),'_BCCAQ2-PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  clim.ncs <- nc_open(clim.file,write=TRUE)
  common.lat <- ncvar_get(clim.ncs,'lat')

##---------------------------------------------------------------------------
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

  pr.input <- paste0('pr_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  print(paste0('Copy ',pr.input))
  file.copy(paste0(data.dir,"/",pr.input),tmp.dir,overwrite=TRUE)

  print('TX opening')
  tasmax.file <- paste0('tasmax_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmax.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
  tasmax.dates <- netcdf.calendar(tasmax.nc)

  print('TN Opening')
  tasmin.file <- paste0('tasmin_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  tasmin.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)
  tasmin.dates <- netcdf.calendar(tasmin.nc)

  print('PR Opening')
  pr.file <- paste0('pr_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  pr.nc <- nc_open(paste0(tmp.dir,pr.file),write=FALSE)
  pr.dates <- netcdf.calendar(pr.nc)
  
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
    pr.subset <- ncvar_get(pr.nc,'pr',start=c(1,i,1),count=c(-1,1,-1))

    flag <- is.na(tasmax.subset[,1])      
    tasmax.list <- vector(mode='list',length=n.lon)
    tasmin.list <- vector(mode='list',length=n.lon)
    pr.list <- vector(mode='list',length=n.lon)

    rtm <- proc.time()     
    tasmax.list <- lapply(seq_len(nrow(tasmax.subset)), function(i) tasmax.subset[i,])
    rm(tasmax.subset)
    tasmin.list <- lapply(seq_len(nrow(tasmin.subset)), function(i) tasmin.subset[i,])     
    rm(tasmin.subset)
    pr.list <- lapply(seq_len(nrow(pr.subset)), function(i) pr.subset[i,])
    rm(pr.subset)
    gc()
    print('Lapply convert to list')
    print(proc.time()-rtm) 

    ##----------------------------------------------------------
    ##Climdex
    rtm <- proc.time()     
      climdex.extremes.for.model(clim.name,clim.ncs,lat.ix,n.lon,flag,
                                 tasmax.list,tasmax.dates,
                                 tasmin.list,tasmin.dates,
                                 pr.list,pr.dates,base=c(1971,2000))
    print('Climdex Time')
    print(proc.time()-rtm) 

    ##----------------------------------------------------------
    rm(tasmax.list)
    rm(tasmin.list)
    rm(pr.list)
    gc()
    print('Lon loop time')
    print(proc.time()-ltm)
  }##Longitude Loop
  nc_close(tasmax.nc)
  nc_close(tasmin.nc)
  nc_close(pr.nc)

  print('Removing lat band files')
  file.remove(paste0(tmp.dir,"/",tasmax.input))
  file.remove(paste0(tmp.dir,"/",tasmin.input))
  file.remove(paste0(tmp.dir,"/",pr.input))
  gc()
}##Latitude File Loop

##Move back
write.dir <- paste0(base.dir,gcm)
move.back <- paste("rsync -av ",tmp.dir,scenario,"/",type," ",write.dir,"/",scenario ,sep='')
##print(move.back)
##system(move.back)

##file.copy(from=paste0(tmp.dir,scenario,"/",type,"/"),to=paste0(write.dir,'/',scenario,'/'),overwrite=TRUE,recursive=TRUE)
file.copy(from=paste0(tmp.dir,scenario,"/",type,"/",template.file),to=paste0(write.dir,'/',scenario,'/',type,'/'),overwrite=TRUE)

nc_close(clim.ncs)


clean.up <- paste("rm ",tmp.dir,scenario,"/",type,"/*",sep='')
print(clean.up)
system(clean.up)

clean.up <- paste("rmdir ",tmp.dir,scenario,"/",type,sep='')
print(clean.up)
system(clean.up)

print("Total Elapsed Time:")
print(proc.time()-ptm)

