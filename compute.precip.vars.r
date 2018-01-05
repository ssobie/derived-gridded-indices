##Script to calculate the degree-day indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file
Rprof('prs.out')
##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

library(ncdf4)
library(PCICt)

library(doParallel)
registerDoParallel(cores=2) 

##--------------------------------
##Degree Day Values
vot11 <- function(pr,thresh) {
  ##prs <- pr[60:300]
  len <- length(pr)
  prs <- pr[c(1:60,300:len)]
  prs[prs < thresh] <- 0
  lens <- rle(prs)
  seps <- lens$lengths[lens$values==0]
  days <- sum(seps[seps <= 5])
}


##Time between R95 Days
prdist <- function(pr,thresh) {
  len <- length(pr)
  prs <- pr[c(1:60,300:len)]
  prs[prs < thresh] <- 0
  lens <- rle(prs)
  seps <- lens$lengths[lens$values==0]
  days <- mean(seps,na.rm=T)
}


##pr.sep <- function(pr,thresh,fac) {
##    yrs <- as.numeric(tapply(pr,fac,function(x,y){vot11(x,y)},thresh))
##}

pr.sep <- function(pr,thresh,fac) {
    yrs <- as.numeric(tapply(pr,fac,function(x,y){prdist(x,y)},thresh))
}


create.base.files <- function(pr.name,gcm,scenario,type=NULL,
                              past.int,proj.int,new.int,
                              data.dir,write.dir) {

  ##files <- list.files(path=data.dir,pattern=gcm,full.name=TRUE)
  ##scen.all <- list.files(path=paste(data.dir,gcm,sep=''),pattern='pr_day',full.name=TRUE)

  files.all <- list.files(path=paste(data.dir,gcm,sep=''),pattern='pr_gcm_prism',full.name=TRUE)

  past.file <- files.all[grep(past.int,files.all)]
  proj.file <- files.all[grep(proj.int,files.all)]

  file.split <- strsplit(past.file,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  write.clim.name <- paste(pr.name,'_annual_',gcm,'_',scenario,'_',run,'_',new.int,'.nc',sep='')
  write.dir <- paste(write.dir,scenario,'/climdex/',gcm,'/',sep='')

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

  var.geog <- ncvar_def(pr.name, units='days', dim=list(x.geog, y.geog, t.geog),
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
  ncatt_put(file.nc,varid=pr.name,attname='units',attval='days')
  ncatt_put(file.nc,varid=pr.name,attname='_FillValue',attval=-32768)
  ncatt_put(file.nc,varid=pr.name,attname='standard_name',attval=toupper(pr.name))
  ncatt_put(file.nc,varid=pr.name,attname='long_name',attval=toupper(pr.name))

  nc_close(file.nc)
  nc_close(nc)
}


##---------------------------------------------------------------

pr.vars.for.model <- function(gcm,scenario,interval,type=NULL,
                              pr.name,store.data,
                              past.int,proj.int,new.int,
                              data.dir,write.dir) {

  pr.files <- list.files(path=paste(data.dir,gcm,sep=''),pattern='pr_gcm_prism',full.name=TRUE)
  pr.past.file <- pr.files[grep(past.int,pr.files)]
  pr.proj.file <- pr.files[grep(proj.int,pr.files)]
 
  file.split <- strsplit(pr.past.file,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  hist.dir <- paste(write.dir,scenario,'/climdex/',gcm,'/',sep='')    

  write.file <- list.files(path=hist.dir,pattern=pr.name,full.name=TRUE)
  print(write.file)
  clim.nc <- nc_open(write.file,write=TRUE)

  ##--------------------------------------------------------------

  pr.past.nc <- nc_open(pr.past.file,write=FALSE)
  pr.proj.nc <- nc_open(pr.proj.file,write=FALSE)


  pr.units <- ncatt_get(pr.past.nc,'pr')$units
  
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
  
  yearly.fac <- as.factor(format(pr.dates,'%Y'))

  ##--------------------------------------------------------------
  ##Compute pr day values and load into newly created netcdf
  for (i in 1:n.lon) {
    print(paste('Lon: ',i,' in ',n.lon,sep=''))
    store.subset <- store.data[i,]

    ##All NA Values
    if (sum(is.na(store.subset)) == length(store.subset)) {
       print('All NA Values')
       pr.list <- vector(mode='list',length=n.lat)
       pr.values <- lapply(pr.list,function(x){return(as.numeric(rep(NA,length(levels(yearly.fac)))))})
    } else { 
  
      pr.past.subset <- ncvar_get(pr.past.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))
      pr.proj.subset <- ncvar_get(pr.proj.nc,'pr',start=c(i,1,1),count=c(1,-1,-1))
      pr.subset <- cbind(pr.past.subset,pr.proj.subset)

##      flag <- is.na(pr.subset[,1])
      pr.list <- list()
      th.list <- list()

      for (j in 1:n.lat) {
        pr.list[[j]] <- pr.subset[j,]
        th.list[[j]] <- store.subset[j]
      }      
    ##No NA Values      
      if (sum(is.na(store.subset)) == 0) {            
           print('No NA Values')
           pr.values <- foreach(
                       pr=pr.list,
                       th=th.list,       
                       .export=c('yearly.fac')
                       ) %dopar% { 
                            pr.values <- pr.sep(pr,th,yearly.fac)
                       }
       
      } else {
         print('Some NA Values')
         ix <- is.na(store.subset)
         pr.values <- vector(mode='list',length=n.lat)
         for (ixf in which(ix)) {
           pr.values[[ixf]] <- as.numeric(rep(NA,length=length(levels(yearly.fac))))
         }
         pr.sub <- pr.list[!ix]
         th.sub <- th.list[!ix]

         pr.subset <- foreach(
                       prs=pr.sub,
                       ths=th.sub,       
                       .export=c('yearly.fac')
                       ) %dopar% { 
                            result <- pr.sep(pr=prs,thresh=ths,fac=yearly.fac)
                       }
         iw <- which(!ix)                       
         for (k in seq_along(iw)) {
           pr.values[[iw[k]]] <- pr.subset[[k]]
         }
      }   
   }
##   pr.values <- c()
##    for (j in 1:n.lat) {
##      print(j)
##      if (is.na(th.list[[j]])) {
##         pr.values[[j]] <- rep(NA,150)
##      } else {
##         pr.values[[j]] <- pr.sep(pr.list[[j]],th.list[[j]],yearly.fac)
##       }
##    }
      
    ncol <- length(pr.values[[1]])         
    pr.matrix <- matrix(unlist(pr.values),nrow=n.lat,ncol=ncol,byrow=TRUE)

##    pr.matrix[flag,] <- NA

    ncvar_put(clim.nc,varid=pr.name,vals=pr.matrix,
              start=c(i,1,1),count=c(1,-1,-1))

  }##Longitude Loop  
  nc_close(clim.nc)
}##Function end

##**************************************************************************************

pr.names <- 'climdex.r95sep'
  
  gcm.list <- c('CanESM2',
                'CCSM4',
                'CNRM-CM5',
                'CSIRO-Mk3-6-0',
                'GFDL-ESM2G',
                'HadGEM2-ES',
                'MIROC5',
                'MPI-ESM-LR',
                'MRI-CGCM3')


##gcm.list <- 'CCSM4'
##gcm.list <- c('MPI-ESM-LR',
##              'MIROC5')
##**************************************************************************************
##-----------------------------------------------------------------------------
##Climdex from the 10km BCCAQ-GCM output
run.bccaq.raw <- function() {
##Running BC only versions
  
  ##data.dir <-  '/home/data/climate/downscale/CMIP5/BCCAQ/'
  data.dir <-  '/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_bc_subset/'
  write.dir <- '/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm/'
  
  scenario <- 'rcp85'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'
  
  pr.names <- sort(pr.names)

  for (gcm in gcm.list) {
    print(gcm)

    first <- lapply(pr.names,create.base.files,
                    gcm,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    data.dir,write.dir)

    second <- pr.days.for.model(gcm,scenario,interval,type='bccaq',
                                    pr.names,
                                    past.int,proj.int,new.int,
                                    data.dir,write.dir)
  }  
}

##-----------------------------------------------------------------------------
##Pr Days from the 800m PRISM adjusted BCCAQ
run.bccaq.prism <- function() {
##Running BC only versions
  
  scenario <- 'rcp85'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  new.int <- '1951-2100'
  region <- 'van_whistler'

  data.dir <- paste('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/',sep='')
  write.dir <- paste('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/',sep='')

  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }

  tmp.base <- tmpdir ##'/local_temp/ssobie/van_whistler/'


  ##Move data to local storage for better I/O
  tmp.dir <- paste('/local_temp/ssobie/',region,'/',sep='')
  if (!file.exists(tmp.dir))
    dir.create(tmp.dir,recursive=TRUE)    
    
  pr.name <- 'r95distETCCDI'
  gcm.list <- c('CSIRO-Mk3-6-0',
                'GFDL-ESM2G',
                'HadGEM2-CC',
                'HadGEM2-ES',
                'inmcm4',
                'MRI-CGCM3')
  gcm.list <- c('MIROC5','MPI-ESM-LR')
                

  print(gcm.list)
  for (gcm in gcm.list) {
    print(gcm)

    move.to <- paste("rsync -av ",data.dir,gcm,"/*pr_gcm_prism* ",tmp.dir,gcm,sep='')
    ##move.to <- paste("rsync -av ",data.dir,gcm,"/pr_day* ",tmp.dir,gcm,sep='')
    print(move.to)
    system(move.to)

    store.file <- paste0(data.dir,'rcp85/climdex/',gcm,'/r95storeETCCDI_ann_BCCAQ-PRISM_',gcm,'_rcp85_r3i1p1_1951-2100.nc')
    store.nc <- nc_open(store.file)
    store.sep <- ncvar_get(store.nc,'r95storeETCCDI',start=c(1,1,1),count=c(-1,-1,1))
    nc_close(store.nc)
    print('Create new files')
    first <- lapply(pr.name,create.base.files,
                    gcm,scenario,type='bccaq',
                    past.int,proj.int,new.int,
                    tmp.dir,tmp.dir)

    print('Calculate Pr Days')
    second <- pr.vars.for.model(gcm,scenario,interval,type='bccaq',
                                pr.name,store.sep,
                                past.int,proj.int,new.int,
                                tmp.dir,tmp.dir)

    move.back <- paste("rsync -av ",tmp.dir,scenario,"/climdex/",gcm," ",write.dir,scenario,"/climdex",sep='')
    print(move.back)
    system(move.back)

    clean.up <- paste("rm ",tmp.dir,gcm,"/*nc" ,sep='')
    print(clean.up)      
    system(clean.up)

    clean.up.dd <- paste("rm ",tmp.dir,scenario,"/climdex/",gcm,"/*nc" ,sep='')
    print(clean.up.dd)      
    system(clean.up.dd)

  }  
}

##**************************************************************************************

run.bccaq.prism()
##run.bccaq.raw()


Rprof(NULL)