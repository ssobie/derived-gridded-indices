##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

rcp26.list <- c('CanESM2',
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

make.netcdf.file <- function(file.name,var.name,clim.data) {

    nc <- nc_open(file.name,write=TRUE)
    ncvar_put(nc,var.name,clim.data)  
    ncatt_put(nc,varid=0,attname='history',attval='')
    ncatt_put(nc,varid=0,attname='driving_experiment',
                      attval='PCIC 12 Ensemble Average, historical+rcp85')
    ncatt_put(nc,varid=0,attname='driving_model_id',
                      attval='PCIC 12 Ensemble Average')
    nc_close(nc)    
}

avg.sub.time <- function(day.file,seas.file,var.name,interval,read.dir,write.dir,fxn) {
  print(paste(read.dir,day.file,sep=''))             
  nc <- nc_open(paste(read.dir,day.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400
  years <- format(time.series,'%Y')
  yrs <- strsplit(interval,'-')[[1]]
  
  st <- head(grep(yrs[1],years),1)
  en <- tail(grep(yrs[2],years),1)
  if (length(en)==0) {
    en <- length(years)
  }   
  len <- en-st+1
  data <- ncvar_get(nc,var.name,start=c(1,1,st),count=c(-1,-1,len))

  yr.fac <- as.factor(years[st:en])
  if (grepl(var.name,'(txxETCCDI_mon|tnxETCCDI_mon|txnETCCDI_mon|tnnETCCDI_mon|rx1dayETCCDI_mon|rx5dayETCCDI_mon)')) {
    clim <- apply(data,c(1,2),function(x,y){mean(tapply(x,y,fxn),na.rm=T)},yr.fac)
  } else {
    clim <- apply(data,c(1,2),fxn,na.rm=TRUE)
  }

  clim[is.infinite(clim)] <- NA
  clim[is.nan(clim)] <- NA
  nc_close(nc)
  return(clim)

}


create.degree.day.climatologies <- function() {
  ds.type <- 'bccaq'
  var.list <- 's30' ##c('cdd','gdd','hdd','fdd') ##'s30' ##

  past.int <- '1971-2000'
  scenario <- 'rcp85'
  region <- 'van_whistler'  
  proj.dir <-  paste0('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/rcp85/degree_days/')
  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/gis_files/degree_days/')
  if (!file.exists(write.dir))
    dir.create(write.dir)

  lonc <- 339
  latc <- 171 

  data.past <- array(NA,c(lonc,latc,12))
##  data.2020s <- array(NA,c(lonc,latc,12))
  data.2050s <- array(NA,c(lonc,latc,12))
##  data.2080s <- array(NA,c(lonc,latc,12))

  for (var.name in var.list) {
    print(var.name)

    ##New files 
    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_annual_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
                                      write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_1971-2000.nc ',sep=''))
##    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_annual_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
##                                      write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_2011-2040.nc ',sep=''))
    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_annual_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
                                      write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_2041-2070.nc ',sep=''))
##    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_annual_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
##                                      write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_2071-2100.nc ',sep=''))

    for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)
      read.dir <- paste(proj.dir,gcm,'/',sep='')
      dd.file <- list.files(path=read.dir,pattern=paste(var.name,'_annual',sep=''))

      data.past[,,g] <- avg.sub.time(dd.file,write.past,var.name,interval=past.int,read.dir,write.dir,mean)
##      data.2020s[,,g] <- avg.sub.time(dd.file,write.proj,var.name,interval='2011-2040',read.dir,write.dir,mean)
      data.2050s[,,g] <- avg.sub.time(dd.file,write.proj,var.name,interval='2041-2070',read.dir,write.dir,mean)
##      data.2080s[,,g] <- avg.sub.time(dd.file,write.proj,var.name,interval='2071-2100',read.dir,write.dir,mean)
    }##GCM Loop   

    clim.past <- apply(data.past,c(1,2),mean,na.rm=TRUE)  
##    clim.2020s <- apply(data.2020s,c(1,2),mean,na.rm=TRUE)  
    clim.2050s <- apply(data.2050s,c(1,2),mean,na.rm=TRUE)  
##    clim.2080s <- apply(data.2080s,c(1,2),mean,na.rm=TRUE)  

    file.past <- paste0(write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_1971-2000.nc')
    make.netcdf.file(file.past,var.name,clim.past)  
##    file.2020s <- paste0(write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_2011-2040.nc')
##    make.netcdf.file(file.2020s,var.name,clim.2020s)  
    file.2050s <- paste0(write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_2041-2070.nc')
    make.netcdf.file(file.2050s,var.name,clim.2050s)  
##    file.2080s <- paste0(write.dir,var.name,'_ensemble_clim_rcp85_r1i1p1_2071-2100.nc')
##    make.netcdf.file(file.2080s,var.name,clim.2080s)  
    
  }##Variable Loop
}

get.clim.fxn <- function(var.name) {

  ann.var.list <- c('fdETCCDI_ann',
                    'suETCCDI_ann',
                    'idETCCDI_ann',
                    'trETCCDI_ann',
                    'gslETCCDI_ann',
                    'wsdiETCCDI_ann',
                    'csdiETCCDI_ann',
                    'sdiiETCCDI_ann',
                    'r10mmETCCDI_ann',
                    'r20mmETCCDI_ann',
                    'cddETCCDI_ann',
                    'cwdETCCDI_ann',
                    'r95pETCCDI_ann',
                    'r99pETCCDI_ann',
                    'prcptotETCCDI_ann',
                    'r95daysETCCDI_ann',
                    'r99daysETCCDI_ann')

  if (any(grepl(var.name,ann.var.list))) {                    
    fxn <- mean
  } else {
    fxn <- switch(var.name,
                  txxETCCDI_mon=max,
                  tnxETCCDI_mon=max,
                  txnETCCDI_mon=min,
                  tnnETCCDI_mon=min,
                  tn10pETCCDI_mon=mean,
                  tx10pETCCDI_mon=mean,
                  tn90pETCCDI_mon=mean,
                  tx90pETCCDI_mon=mean,
                  dtrETCCDI_mon=mean,
                  rx1dayETCCDI_mon=max,
                  rx5dayETCCDI_mon=max)

  }

  return(fxn)

}


create.climdex.climatologies <- function() {
  ds.type <- 'bccaq'
  var.list <- c('fdETCCDI_ann',
                    'suETCCDI_ann',
                    'idETCCDI_ann',
                    'trETCCDI_ann',
                    'gslETCCDI_ann',
                    'wsdiETCCDI_ann',
                    'csdiETCCDI_ann',
                    'sdiiETCCDI_ann',
                    'r10mmETCCDI_ann',
                    'r20mmETCCDI_ann',
                    'cddETCCDI_ann',
                    'cwdETCCDI_ann',
                    'r95pETCCDI_ann',
                    'r99pETCCDI_ann',
                    'prcptotETCCDI_ann',
                    'r95daysETCCDI_ann',
                    'r99daysETCCDI_ann',
                    'tn10pETCCDI_mon',
                    'tx10pETCCDI_mon',
                    'tn90pETCCDI_mon',
                    'tx90pETCCDI_mon',
                    'dtrETCCDI_mon')


 var.list <- c('txxETCCDI_mon',
                    'tnxETCCDI_mon',
                    'txnETCCDI_mon',
                    'tnnETCCDI_mon',
                    'rx1dayETCCDI_mon',
                    'rx5dayETCCDI_mon')
  var.names <- unlist(lapply(strsplit(var.list,'_'),function(x){return(x[1])}))


  past.int <- '1971-2000'
  scenario <- 'rcp85'
  
  ##proj.dir <-  '/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_south_island_subset/rcp85/climdex/'
  ##write.dir  <- '/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_south_island_subset/gis_files/climdex/'

region <- 'nanaimo'

  proj.dir <-  paste0('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/rcp85/climdex/')
  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/gis_files/climdex/')
  type <- 'BCCAQ-PRISM'
  if (!file.exists(write.dir))
    dir.create(write.dir)

  lonc <- 73
  latc <- 58 

  data.past <- array(NA,c(lonc,latc,12))
  data.2020s <- array(NA,c(lonc,latc,12))
  data.2050s <- array(NA,c(lonc,latc,12))
  data.2080s <- array(NA,c(lonc,latc,12))

  for (v in seq_along(var.list)) {
    var.title <- var.list[v]
    var.name <- var.names[v]
    print(var.name)

    ##New files 
    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ-PRISM_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_1971-2000.nc ',sep=''))
    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ-PRISM_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2011-2040.nc ',sep=''))
    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ-PRISM_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2041-2070.nc ',sep=''))
    system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ-PRISM_ACCESS1-0_rcp85_r1i1p1_1951-2100.nc ',
                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2071-2100.nc ',sep=''))

    for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)
      read.dir <- paste(proj.dir,gcm,'/',sep='')
      dd.file <- list.files(path=read.dir,pattern=paste(var.title,'_',sep=''))
      
      data.past[,,g] <- avg.sub.time(dd.file,write.past,var.name,interval=past.int,read.dir,write.dir,get.clim.fxn(var.title))
      data.2020s[,,g] <- avg.sub.time(dd.file,write.proj,var.name,interval='2011-2040',read.dir,write.dir,get.clim.fxn(var.title))
      data.2050s[,,g] <- avg.sub.time(dd.file,write.proj,var.name,interval='2041-2070',read.dir,write.dir,get.clim.fxn(var.title))
      data.2080s[,,g] <- avg.sub.time(dd.file,write.proj,var.name,interval='2071-2100',read.dir,write.dir,get.clim.fxn(var.title))
    }##GCM Loop   
    
    clim.past <- apply(data.past,c(1,2),mean,na.rm=TRUE)  
    clim.2020s <- apply(data.2020s,c(1,2),mean,na.rm=TRUE)  
    clim.2050s <- apply(data.2050s,c(1,2),mean,na.rm=TRUE)  
    clim.2080s <- apply(data.2080s,c(1,2),mean,na.rm=TRUE)  

    file.past <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_1971-2000.nc')
    make.netcdf.file(file.past,var.name,clim.past)  
    file.2020s <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2011-2040.nc')
    make.netcdf.file(file.2020s,var.name,clim.2020s)  
    file.2050s <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2041-2070.nc')
    make.netcdf.file(file.2050s,var.name,clim.2050s)  
    file.2080s <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2071-2100.nc')
    make.netcdf.file(file.2080s,var.name,clim.2080s)  
    
  }##Variable Loop
}

create.return.period.climatologies <- function() {
  ds.type <- 'bccaq'
##  var.list <- c('pr_RPCI20','tasmax_RPCI20','tasmin_RPCI20')
##  var.names <- c('pr_RP20','tasmax_RP20','tasmin_RP20')
  var.list <- c('rx2dayETCCDI_RPCI10','rx5dayETCCDI_RPCI10')
  var.names <- c('rx2day_RP10','rx5day_RP10')

  past.int <- '1971-2000'
  scenario <- 'rcp85'
  region <- 'van_whistler'
  proj.dir <-  paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/rcp85/return_periods/')
  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_',region,'_subset/gis_files/return_periods/')
  if (!file.exists(write.dir))
    dir.create(write.dir)

  lonc <- 339
  latc <- 171 

  data.past <- array(NA,c(lonc,latc,12))
##  data.2020s <- array(NA,c(lonc,latc,12))
  data.2050s <- array(NA,c(lonc,latc,12))
##  data.2080s <- array(NA,c(lonc,latc,12))

  for (v in seq_along(var.list)) {
    var.title <- var.list[v]
    var.name <- var.names[v]
    print(var.name)

    ##New files 
    system(paste('cdo -s timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ_PRISM_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_1971-2000.nc ',sep=''))
##    system(paste('cdo -s timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ_PRISM_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
##                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2011-2040.nc ',sep=''))
    system(paste('cdo -s timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ_PRISM_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2041-2070.nc ',sep=''))
##    system(paste('cdo -s timmean ',proj.dir,'ACCESS1-0/',var.title,'_BCCAQ_PRISM_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
##                                      write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2071-2100.nc ',sep=''))

    for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)
      read.dir <- paste(proj.dir,gcm,'/',sep='')
      rp.files <- list.files(path=read.dir,pattern=paste(var.title,'_BCCAQ_PRISM',sep=''),full.name=TRUE)
      past.file <- rp.files[grep('1971-2000',rp.files)]
      nc.past <- nc_open(past.file)
      data.past[,,g] <- ncvar_get(nc.past,'rp.10')[,,2]
##      file.2020s <- rp.files[grep('2011-2040',rp.files)]
##      nc.2020s <- nc_open(file.2020s)
##      data.2020s[,,g] <- ncvar_get(nc.2020s,'rp.10')[,,2]
      file.2050s <- rp.files[grep('2041-2070',rp.files)]
      nc.2050s <- nc_open(file.2050s)
      data.2050s[,,g] <- ncvar_get(nc.2050s,'rp.10')[,,2]
##      file.2080s <- rp.files[grep('2071-2100',rp.files)]
##      nc.2080s <- nc_open(file.2080s)
##      data.2080s[,,g] <- ncvar_get(nc.2080s,'rp.10')[,,2] 

    }##GCM Loop   

    clim.past <- apply(data.past,c(1,2),mean,na.rm=TRUE)  
##    clim.2020s <- apply(data.2020s,c(1,2),mean,na.rm=TRUE)  
    clim.2050s <- apply(data.2050s,c(1,2),mean,na.rm=TRUE)  
##    clim.2080s <- apply(data.2080s,c(1,2),mean,na.rm=TRUE)  

    file.past <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_1971-2000.nc')
    make.netcdf.file(file.past,'rp.10',clim.past)  
##    file.2020s <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2011-2040.nc')
##    make.netcdf.file(file.2020s,'rp.20',clim.2020s)  
    file.2050s <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2041-2070.nc')
    make.netcdf.file(file.2050s,'rp.10',clim.2050s)  
##    file.2080s <- paste0(write.dir,var.name,'_ensemble_',region,'_rcp85_r1i1p1_2071-2100.nc')
##    make.netcdf.file(file.2080s,'rp.20',clim.2080s)  
    
  }##Variable Loop
}

 
create.seasonal.climatologies <- function() {
  ds.type <- 'bccaq'
  var.list <- c('pr','tasmax','tasmin')

  type.list <- c('winter','spring','summer','fall') ##'annual' ##
  base.type <- 'seas' ##'ann' ##
  region <- 'bc'
  proj.dir <-  paste0('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/rcp85/seasonal/')
  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2/high_res_downscaling/bccaq_gcm_',region,'_subset/gis_files/seasonal/')
  past.int <- '1971-2000'
  scenario <- 'rcp85'  

  if (!file.exists(write.dir))
    dir.create(write.dir)

  lonc <- 313
  latc <- 145 

  data.past <- array(NA,c(lonc,latc,12))
##  data.2020s <- array(NA,c(lonc,latc,12))
  data.2050s <- array(NA,c(lonc,latc,12))
##  data.2080s <- array(NA,c(lonc,latc,12))

  for (s in seq_along(type.list)) {
    type <- type.list[s]
    for (var.name in var.list) {
      print(var.name)

      ##New files 
      if (base.type=='ann') {
        system(paste('cp ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_1971-2000.nc ',sep=''))
##        system(paste('cp ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
##                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2011-2040.nc ',sep=''))
        system(paste('cp ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2041-2070.nc ',sep=''))
##        system(paste('cp ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
##                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2071-2100.nc ',sep=''))
      } else {
        system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_1971-2000.nc ',sep=''))
##        system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
##                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2011-2040.nc ',sep=''))
        system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2041-2070.nc ',sep=''))
##        system(paste('cdo -s -O timmean ',proj.dir,'ACCESS1-0/',var.name,'_',base.type,'_BCCAQ_ACCESS1-0_rcp85_r1i1p1_1971-2000.nc ',
##                                          write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2071-2100.nc ',sep=''))
      }

      for (g in seq_along(gcm.list)) {
        gcm <- gcm.list[g]
        print(gcm)
        read.dir <- paste(proj.dir,gcm,'/',sep='')
        type.files <- list.files(path=read.dir,pattern=paste(var.name,'_',base.type,sep=''),full.name=TRUE)
        past.file <- type.files[grep('1971-2000',type.files)]
        nc.past <- nc_open(past.file)

        data.past[,,g] <- ncvar_get(nc.past,var.name)[,,s]
##        file.2020s <- type.files[grep('2011-2040',type.files)]
##        nc.2020s <- nc_open(file.2020s)
##        data.2020s[,,g] <- ncvar_get(nc.2020s,var.name)[,,s]
        file.2050s <- type.files[grep('2041-2070',type.files)]
        nc.2050s <- nc_open(file.2050s)
        data.2050s[,,g] <- ncvar_get(nc.2050s,var.name)[,,s]
##        file.2080s <- type.files[grep('2071-2100',type.files)]
##        nc.2080s <- nc_open(file.2080s)
##        data.2080s[,,g] <- ncvar_get(nc.2080s,var.name)[,,s] 
      }##GCM Loop   

      clim.past <- apply(data.past,c(1,2),mean,na.rm=TRUE)  
##      clim.2020s <- apply(data.2020s,c(1,2),mean,na.rm=TRUE)  
      clim.2050s <- apply(data.2050s,c(1,2),mean,na.rm=TRUE)  
##      clim.2080s <- apply(data.2080s,c(1,2),mean,na.rm=TRUE)  

      file.past <- paste0(write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_1971-2000.nc')
      make.netcdf.file(file.past,var.name,clim.past)  
##      file.2020s <- paste0(write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2011-2040.nc')
##      make.netcdf.file(file.2020s,var.name,clim.2020s)  
      file.2050s <- paste0(write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2041-2070.nc')
      make.netcdf.file(file.2050s,var.name,clim.2050s)  
##      file.2080s <- paste0(write.dir,var.name,'_ensemble_',type,'_',region,'_rcp85_r1i1p1_2071-2100.nc')
##      make.netcdf.file(file.2080s,var.name,clim.2080s)  
    
    }##Variable Loop
  }##Season Loop
}
