##Script to calculate and write the standard set of derived variables
##for the 800m data

##Rprof('cedar.derived.out')

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

new.tasmin.climdex.extremes.for.model <- function(climdex.names,clim.ncs,lat.ix,n.lon,flag,base,
                                              tasmin.list,tasmin.dates) {

   climdex.matrix <- matrix(NA,nrow=n.lon,ncol=150)
   climdex.list <- rep(list(climdex.matrix),length(climdex.names))
   rm(climdex.matrix)

   flen <- sum(!flag)
   tasmin.sub <- tasmin.list[!flag]
   rm(tasmin.list)
   sub.indices  <- which(!flag)

   sqnce <- seq(0,flen,by=200)
   sx <- findInterval(flen,sqnce)
   itvls <- c(sqnce[1:sx],flen-sqnce[sx]+sqnce[sx])

   for (s in 1:sx) {
      st <- itvls[s] + 1
      en <- itvls[s+1]
      len <- en-st+1
      tasmin.slice <- tasmin.sub[st:en]
      sub.ix <- sub.indices[st:en]
      climdex.objects <- foreach(
                                 tasmin=tasmin.slice,
                                 .export=c('climdexInput.raw','tasmin.dates','base')
                                 ) %dopar% {
                                   objects <- climdexInput.raw(tmin=tasmin,
                                                               tmin.dates=tasmin.dates,
                                                               base=base)
                                 }
      rm(tasmin.slice)                                       
                                         
      for (k in seq_along(climdex.names)) {
        clim.name <- climdex.names[k]        
        clim.nc <- clim.ncs[[k]]
        clim.fx <- clim.fxns[[clim.name]]
        climdex.values <- foreach(
                                 obj=climdex.objects,
                                 .export='clim.fx'
                                 ) %dopar% {
                                 climdex.values <- clim.fx(obj)
                                 }
        ncol <- length(climdex.values[[1]])
        sub.matrix <- matrix(unlist(climdex.values),nrow=len,ncol=ncol,byrow=TRUE)       
        rm(climdex.values)
        climdex.list[[k]][sub.ix,] <- sub.matrix
        rm(sub.matrix)
      }##Climdex names vars loop
      rm (climdex.objects)
      gc()
   }

   for (k in seq_along(climdex.names)) {
        clim.nc <- clim.ncs[[k]]
        clim.name <- climdex.names[k]        
        climdex.info <- get.climdex.info(clim.name)
        ncvar_put(clim.nc,varid=climdex.info[1],vals=climdex.list[[k]],
                  start=c(1,lat.ix,1),count=c(-1,1,-1))
   }
   rm(climdex.list)
   gc()
}


calc.climdex.vars <- function(climdex.objects,climdex.names,climdex.ncs,n.lon,flen,flag) {

   for (k in seq_along(climdex.names)) {
        clim.name <- climdex.names[k]        
        clim.nc <- clim.ncs[[k]]

        clim.fx <- clim.fxns[[clim.name]]
        climdex.info <- get.climdex.info(clim.name)

        climdex.values <- foreach(
                                 obj=climdex.objects,
                                 .export='clim.fx'
                                 ) %dopar% {
                                 climdex.values <- clim.fx(obj)
                                 }
        ncol <- length(climdex.values[[1]])
        climdex.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
        sub.matrix <- matrix(unlist(climdex.values),nrow=flen,ncol=ncol,byrow=TRUE)
        rm(climdex.values)
        climdex.matrix[!flag,] <- sub.matrix
        rm(sub.matrix)
        ncvar_put(clim.nc,varid=climdex.info[1],vals=climdex.matrix,
                  start=c(1,lat.ix,1),count=c(-1,1,-1))
        rm(climdex.matrix)                  
    }##Climdex names vars loop

}


tasmax.climdex.extremes.for.model <- function(climdex.names,clim.ncs,lat.ix,n.lon,flag,base,
                                              tasmax.list,tasmax.dates) {

   flen <- sum(!flag)
   tasmax.sub <- tasmax.list[!flag]
   
   climdex.objects <- foreach(
                             tasmax=tasmax.sub,
                             .export=c('climdexInput.raw','tasmax.dates','base')
                             ) %dopar% {
                                objects <- climdexInput.raw(tmax=tasmax,
                                                            tmax.dates=tasmax.dates,
                                                            base=base)
                             }
   calc.climdex.vars(climdex.objects,climdex.names,climdex.ncs,n.lon,flen,flag) 
}

tasmin.climdex.extremes.for.model <- function(climdex.names,clim.ncs,lat.ix,n.lon,flag,base,
                                              tasmin.list,tasmin.dates) {

   flen <- sum(!flag)
   tasmin.sub <- tasmin.list[!flag]

   climdex.objects <- foreach(
                             tasmin=tasmin.sub,
                             .export=c('climdexInput.raw','tasmin.dates','base')
                             ) %dopar% {
                                objects <- climdexInput.raw(tmin=tasmin,
                                                            tmin.dates=tasmin.dates,
                                                            base=base)
                             }
   calc.climdex.vars(climdex.objects,climdex.names,climdex.ncs,n.lon,flen,flag) 
   
}


tas.climdex.extremes.for.model <- function(climdex.names,clim.ncs,lat.ix,n.lon,flag,base,
                                              tasmax.list,tasmax.dates,tasmin.list,tasmin.dates) {

   flen <- sum(!flag)
   tasmax.sub <- tasmax.list[!flag]
   tasmin.sub <- tasmin.list[!flag]

   climdex.objects <- foreach(
                             tasmax=tasmax.sub,
                             tasmin=tasmin.sub,
                             .export=c('climdexInput.raw','tasmax.dates','tasmin.dates','base')
                             ) %dopar% {
                                objects <- climdexInput.raw(tmax=tasmax,tmax.dates=tasmax.dates,
                                                            tmin=tasmin,tmin.dates=tasmin.dates,    
                                                            base=base)
                             }
   calc.climdex.vars(climdex.objects,climdex.names,climdex.ncs,n.lon,flen,flag) 
   
}


pr.climdex.extremes.for.model <- function(climdex.names,clim.ncs,lat.ix,n.lon,flag,base,
                                          pr.list,pr.dates) {

   flen <- sum(!flag)
   pr.sub <- pr.list[!flag]

   climdex.objects <- foreach(
                             pr=pr.sub,
                             .export=c('climdexInput.raw','pr.dates','base')
                             ) %dopar% {
                                objects <- climdexInput.raw(tmin=pr,
                                                            tmin.dates=pr.dates,
                                                            base=base)
                             }
   calc.climdex.vars(climdex.objects,climdex.names,climdex.ncs,n.lon,flen,flag) 
   
}



##--------------------------------------------------------------
##****************************************************************

if (1==0) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

tmp.dir <- '/local_temp/ssobie/prism/' ##tmpdir

gcm <- 'CNRM-CM5'
scenario <- 'rcp85'
run <- 'r1i1p1'
interval <- '1951-2100'
type <- 'annual'
varname <- 'tasmin'

base.period <- c(1971,2000)

##---------------------------------------------------------------------------
##Climdex Files
if (type=='annual') {
   if (varname=='tasmax') {
     clim.names <- c('id','su','wsdi')
   }  
   if (varname=='tasmin') {
     clim.names <- c('csdi','fd','tr')
   }
   if (varname=='tas') {
     clim.names <- 'gsl'
   }
   if (varname=='prdays') {
     clim.names <- c('cdd','cwd','prcptot','sdii')
     varmame <- 'pr'
   }
   if (varname=='pr95') {
     clim.names <- c('r95days','r95ptot','r95store')
     varname <- 'pr'
   }        
   if (varname=='pr99') {              
     clim.names <- c('r99days','r99ptot','r99store')
     varname <- 'pr'
   }
} else {
   clim.names <- type
}

climdex.names <- paste0('climdex.',clim.names)

##Latitude Bands
##lat.st <- seq(48.1,59.9,0.1) ##format(seq(48.0,59.9,0.1),nsmall=1)
##lat.en <- seq(48.2,60.0,0.1) ##format(seq(48.1,60.0,0.1),nsmall=1)

lat.st <- seq(48.9,59.9,0.1)
lat.en <- seq(49.0,60.0,0.1)

len <- length(lat.st)

if (cedar) {
   base.dir <- '/scratch/ssobie/prism/'
} else {
   base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
}

data.dir <- paste0(base.dir,gcm,'/lat_split/')
template.dir <- paste0(base.dir,gcm,'/template/',scenario,'/climdex/')
all.files <- list.files(path=template.dir,pattern='ETCCDI')
all.names <- unlist(lapply(strsplit(all.files,'_'),function(x){x[1]}))
template.match <- all.names %in% paste0(clim.names,'ETCCDI')
template.files <- paste0(template.dir,all.files[template.match])

##Move data to local storage for better I/O
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}
dir.create(paste0(tmp.dir,scenario),recursive=T)

rtm <- proc.time()
print('Copying')
print(template.files)

##file.copy(from=template.files,to=paste0(tmp.dir,scenario,'/climdex/'),overwrite=TRUE,recursive=TRUE) ##

print('Move template time')
print(proc.time()-rtm)

  climdex.dir <- paste0(tmp.dir,'/',scenario,'/climdex/')
  clim.ncs <- vector(mode='list',length=length(climdex.names))
  for (d in seq_along(climdex.names)) {
    clim.name <- climdex.names[d]
    climdex.info <- get.climdex.info(clim.name)
    climdex.var <- climdex.info[1]
    climdex.calendar <- climdex.info[2]
    clim.file <- paste0(climdex.dir,climdex.var,'_',tolower(climdex.calendar),'_BCCAQ2-PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
    clim.ncs[[d]] <- nc_open(clim.file,write=TRUE)
  }
  common.lat <- ncvar_get(clim.ncs[[1]],'lat')

##---------------------------------------------------------------------------
##---------------------------------------------------------------------------

##Iterate over the latitude files
for (i in 1:len) {
  ftm <- proc.time()  
  lat.interval <- paste0(sprintf('%.1f',lat.st[i]),'-',sprintf('%.1f',lat.en[i]))

  data.input <- paste0(varname,'_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  print(paste0('Copy ',data.input))
  ##file.copy(paste0(data.dir,"/",data.input),tmp.dir,overwrite=TRUE)

  print('Data opening')
  input.file <- paste0(varname,'_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  input.nc <- nc_open(paste0(tmp.dir,input.file),write=FALSE)
  input.dates <- netcdf.calendar(input.nc)
  
  lon <- ncvar_get(input.nc,'lon')
  lat <- ncvar_get(input.nc,'lat')
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
    input.subset <- ncvar_get(input.nc,varname,start=c(1,i,1),count=c(-1,1,-1))
    flag <- is.na(input.subset[,1])
    input.list <- vector(mode='list',length=n.lon)

    rtm <- proc.time()
    input.list <- lapply(seq_len(nrow(input.subset)), function(i) input.subset[i,])
    rm(input.subset)
    print('Lapply convert to list')
    print(proc.time()-rtm)

    ##Climdex
    if (varname=='tasmax') {
      tasmax.climdex.extremes.for.model(climdex.names,clim.ncs,lat.ix,n.lon,flag,base=base.period,
                                        input.list,input.dates)
    }                                        
    if (varname=='tasmin') {
      new.tasmin.climdex.extremes.for.model(climdex.names,clim.ncs,lat.ix,n.lon,flag,base=base.period,
                                        input.list,input.dates)
    }                                        
    if (grepl('pr',varname)) {
      pr.climdex.extremes.for.model(climdex.names,clim.ncs,lat.ix,n.lon,flag,base=base.period,
                                    input.list,input.dates)
    }                                        

    ##----------------------------------------------------------
    rm(input.list)

    print('Lon loop time')
    print(proc.time()-ltm)
  }##Longitude Loop
  nc_close(input.nc)

  print('Removing lat band files')
  file.remove(paste0(tmp.dir,"/",data.input))

  print('File loop time (time per file)')
  print(proc.time()-ftm)
  browser()
}##Latitude File Loop

##Move back
write.dir <- paste0(base.dir,gcm)
move.back <- paste("rsync -av ",tmp.dir,scenario,"/",type," ",write.dir,"/",scenario ,sep='')
##print(move.back)
##system(move.back)

file.copy(from=paste0(tmp.dir,scenario,"/",type,"/"),to=paste0(write.dir,'/',scenario,'/'),overwrite=TRUE,recursive=TRUE)




  for (d in seq_along(clim.names)) {
    nc_close(clim.ncs[[d]])
  }


##Rprof(NULL)


clean.up <- paste("rm ",tmp.dir,scenario,"/",type,"/*",sep='')
print(clean.up)
##system(clean.up)

clean.up <- paste("rmdir ",tmp.dir,scenario,"/",type,sep='')
print(clean.up)
##system(clean.up)

print("Total Elapsed Time:")
print(proc.time()-ptm)

