##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)
library(raster)

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

create.ensemble.annual.climatologies <- function(var.name,type,season,interval,gcm.list,region) {

  scenario <- 'rcp85'

  ##proj.dir <-   '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
  proj.dir <-   paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/')
  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/',scenario,'/',type,'/ENSEMBLE/')

  print(write.dir)

  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }
  gx <- length(gcm.list)
  lonc <- 3121
  latc <- 1441

  data.past <- array(NA,c(lonc,latc,gx))
  print(var.name)
  ###copy.dir <- paste0(proj.dir,'ACCESS1-0/rcp85/',type,'/climatologies/') ### For BC
    copy.dir <- paste0(proj.dir,'rcp85/',type,'/ACCESS1-0/')
    ##New files     
    copy.files <- list.files(path=copy.dir,pattern=interval)

    var.files <- copy.files[grep(paste0('^',var.name),copy.files)]
    copy.file <- var.files[grep(paste0(season,'_'),var.files)]
    new.file <- gsub('ACCESS1-0','ENSEMBLE',copy.file)

    file.copy(from=paste0(copy.dir,copy.file),
              to=paste0(write.dir,new.file),overwrite=T)
    file.split <- strsplit(copy.file,'_')[[1]]

    data.ens <- c()

    for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      read.dir <- paste(proj.dir,'/rcp85/',type,'/',gcm,'/',sep='')
      all.files <- list.files(path=read.dir,pattern=interval)
      var.files <- all.files[grep(paste0('^',var.name),all.files)]
      seas.file <- var.files[grep(paste0(season,'_'),var.files)]

      print('Selected file')
      print(seas.file)
      box.read <- brick(paste0(read.dir,seas.file))
      
      if (is.null(data.ens)) {
        data.ens <- box.read
      } else {
        data.ens <- stack(data.ens,box.read)
      }
    }##GCM Loop   

    if (grepl('RP',season)) {
      flag <- data.ens==1111
      data.ens[flag] <- NA
    }
    data.tmp <- aperm(as.matrix(calc(data.ens,mean,na.rm=T)),c(2,1))
    data.write <- data.tmp[,ncol(data.tmp):1]

    file.write <- paste0(write.dir,new.file)
    if (season=='RP20') {       
      var.name <- 'rp.20'
    }
    if (season=='RP5') {       
      var.name <- 'rp.5'
    }
    if (season=='RP50') {       
      var.name <- 'rp.50'
    }


    make.netcdf.file(file.write,var.name,data.write)  

}


create.ensemble.seas.mon.climatologies <- function(var.name,type,season,interval,gcm.list,region) {

  scenario <- 'rcp85'

  ###proj.dir <-   '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/' ##For BC
  proj.dir <-   paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/')
  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/',scenario,'/',type,'/ENSEMBLE/')
  print(write.dir)

  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  seas.ix <- switch(season,
                    winter=1,spring=2,summer=3,fall=4,
                    january=1,february=2,march=3,april=4,may=5,june=6,july=7,august=8,september=9,october=10,november=11,december=12)
  gx <- length(gcm.list)
  lonc <- 3121
  latc <- 1441

  data.past <- array(NA,c(lonc,latc,gx))
  print(var.name)
  ###copy.dir <- paste0(proj.dir,'ACCESS1-0/rcp85/',type,'/climatologies/') ##For BC
  copy.dir <- paste0(proj.dir,'rcp85/',type,'/ACCESS1-0/')
    ##New files 
    copy.files <- list.files(path=copy.dir,pattern=interval)
    var.files <- copy.files[grep(var.name,copy.files)]

    if (grepl('(seas|mon)',type)) {
      copy.file <- var.files[grep(type,var.files)]
      inter.file <- gsub('ACCESS1-0','ENSEMBLE',copy.file)
      new.file <- gsub(type,season,inter.file)      
    } else {
      copy.file <- var.files[grep('seasonal',var.files)] ##For seasonal climdex files
      inter.file <- gsub('ACCESS1-0','ENSEMBLE',copy.file)
      new.file <- gsub('seasonal',season,inter.file)      
    }
    file.copy(from=paste0(copy.dir,copy.file),
              to=paste0(write.dir,'tmp.nc'),overwrite=T)
    system(paste0('cdo -s -O timmean ',write.dir,'tmp.nc ',
                                      write.dir,new.file))

    file.split <- strsplit(copy.file,'_')[[1]]
    run <- file.split[grep('r*i1p1',file.split)]

    data.ens <- c()

    for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)
      ###read.dir <- paste(proj.dir,'/',gcm,'/rcp85/',type,'/climatologies/',sep='') ##For BC
      read.dir <- paste(proj.dir,'/rcp85/',type,'/',gcm,'/',sep='')
      all.files <- list.files(path=read.dir,pattern=interval)
      var.files <- all.files[grep(var.name,all.files)]
      if (grepl('(seas|mon)',type)) {
        seas.file <- var.files[grep(type,var.files)]
      } else { ##For climdex files
        seas.file  <- var.files[grep('seasonal',var.files)]
      }

      box.read <- brick(paste0(read.dir,seas.file))
      box.sub <- subset(box.read,seas.ix)

      if (is.null(data.ens)) {
        data.ens <- box.sub
      } else {
        data.ens <- stack(data.ens,box.sub)
      }
    }##GCM Loop   

    data.tmp <- aperm(as.matrix(calc(data.ens,mean)),c(2,1))
    data.write <- data.tmp[,ncol(data.tmp):1]

    file.write <- paste0(write.dir,new.file)
    make.netcdf.file(file.write,var.name,data.write)  

}


create.anomalies.from.climatologies <- function(var.name,type,season,region,proj.int,
                                                past.int='1971-2000') {
  
  read.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/rcp85/',type,'/ENSEMBLE/')
  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/rcp85/',type,'/ENSEMBLE/anomalies/')

  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  var.files <- list.files(path=read.dir,pattern=var.name)
  seas.files <- var.files[grep(season,var.files)]
  past.file <- seas.files[grep(past.int,seas.files)]
  proj.file <- seas.files[grep(proj.int,seas.files)]

  dates.file <- gsub(pattern=past.int,replacement=paste0(past.int,'_',proj.int),past.file)
  anoms.file <- gsub(pattern='climatology',replacement='anomaly',dates.file)
  prct.file <- gsub(pattern='climatology',replacement='percent',dates.file)

  ##Absolute anomalies
  abs.anoms <- paste0('cdo -s -O sub ',read.dir,proj.file,' ',
                                        read.dir,past.file,' ',
                                        write.dir,anoms.file)
  system(abs.anoms)

  do.prct <- grepl('(rx|pr|r1|r2|r9)',var.name)
  print(do.prct)
  if (do.prct) {
    prct.mid <- paste0('cdo -s -O div ',write.dir,anoms.file,' ',
                                        read.dir,past.file,' ',
                                        write.dir,'tmp.nc')
    system(prct.mid)                                        
    prct.anoms <- paste0('cdo -s -O mulc,100 ',write.dir,'tmp.nc ',
                                           write.dir,prct.file)
    system(prct.anoms)
    file.remove(paste0(write.dir,'tmp.nc'))
  }
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

 
if (1==0) {
intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
var.list <- c('pr') ##,'tasmax','tasmin')
seasons <- 'RP50'
type <- 'return_periods'

regions <- c('van_coastal_health') ##,'bella_health','northeast','willow_road')

for (region in regions) {
  for (interval in intervals) {
    for (var.name in var.list) {
      for (season in seasons) {
        create.ensemble.annual.climatologies(var.name,type,season,interval,gcm.list,region) 
      }
    } 
  }
}
}

if (1==0) {
intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
var.list <- c('pr','tasmax','tasmin')
##seasons <- c('winter','spring','summer','fall')
seasons <- c('january','february','march','april','may','june','july','august','september','october','november','december')
##regions <- c('willow_road','van_coastal_health','bella_health','northeast')
regions <- 'bc'

for (region in regions) {
  for (interval in intervals) {
    for (var.name in var.list) {
      for (season in seasons) {
        create.ensemble.seas.mon.climatologies(var.name,type='monthly',season,interval,gcm.list,region) 
      }
    }
  }
}
}


if (1==0) {
intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
var.list <-  c('txxETCCDI','txnETCCDI','tnnETCCDI','tnxETCCDI','dtrETCCDI',
                 'rx1dayETCCDI','rx5dayETCCDI')
               

seasons <- c('winter','spring','summer','fall')
regions <- c('willow_road','van_coastal_health','bella_health','northeast')

for (region in regions) {
  for (interval in intervals) {
    for (var.name in var.list) {
      for (season in seasons) {
        create.ensemble.seas.mon.climatologies(var.name,type='climdex',season,interval,gcm.list,region) 
      }
    }
  }
}
}

if (1==0) {
intervals <- c('2011-2040','2041-2070','2071-2100')
  intervals <- '1971-2000'        
var.list <-  c('fdETCCDI','suETCCDI','su30ETCCDI','idETCCDI','trETCCDI','gslETCCDI',
                'txxETCCDI','txnETCCDI','tnnETCCDI','tnxETCCDI','dtrETCCDI',
                'rx1dayETCCDI','rx5dayETCCDI',
                'sdiiETCCDI','r10mmETCCDI','r20mmETCCDI','cwdETCCDI',
                'cwdETCCDI','cddETCCDI','prcptotETCCDI',
                'r95pETCCDI','r99pETCCDI','r95daysETCCDI','r99daysETCCDI')
seasons <- 'annual'
type <- 'climdex'
options(warn=2)
regions <- c('northeast','bella_health') ##,'van_coastal_health','bella_health','northeast')

for (region in regions) {
  for (var.name in var.list) {
    for (proj.int in intervals) {
      for (season in seasons) {
         create.ensemble.annual.climatologies(var.name,type,season,proj.int,gcm.list,region)                                              
      }
    } 
  }
}
}



if (1==0) {
intervals <- c('2011-2040','2041-2070','2071-2100')
##var.list <-  c('fdETCCDI','suETCCDI','su30ETCCDI','idETCCDI','trETCCDI','gslETCCDI',
##                 'txxETCCDI','txnETCCDI','tnnETCCDI','tnxETCCDI','dtrETCCDI',
##                 'rx1dayETCCDI','rx5dayETCCDI',
##                 'sdiiETCCDI','r10mmETCCDI','r20mmETCCDI','cwdETCCDI','cddETCCDI','prcptotETCCDI',
##                 'r95pETCCDI','r99pETCCDI','r95daysETCCDI','r99daysETCCDI')
var.list <- c('tasmax','tasmax','tasmax',
              'tasmin','tasmin','tasmin')
seas.list <- c('annual_quantile_975','annual_quantile_990','annual_quantile_996',
               'annual_quantile_004','annual_quantile_010','annual_quantile_025')
type <- 'annual_quantiles'

regions <- c('van_coastal_health','bella_health')
for (region in regions) {
    for (v in seq_along(var.list)) {
      var.name <- var.list[v]
      season <- seas.list[v]
      for (proj.int in intervals) {
         create.anomalies.from.climatologies(var.name,type,season,region,proj.int,
                                             past.int='1971-2000')
      } 
    }
}
}

if (1==1) {
intervals <- c('2011-2040','2041-2070','2071-2100')
var.list <- c('r95pETCCDI','r99pETCCDI','r95daysETCCDI','r99daysETCCDI')
seas.list <- 'annual'
type <- 'climdex'
regions <- c('van_coastal_health','bella_health')
for (region in regions) {
    for (var.name in var.list) {
      for (proj.int in intervals) {
         for (season in seas.list) {
           create.anomalies.from.climatologies(var.name,type,season,region,proj.int,
                                               past.int='1971-2000')
         }
      } 
    }
}
}


##Annual Quantiles
if (1==0) {
intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
var.list <- c('tasmax','tasmax','tasmax',
              'tasmin','tasmin','tasmin')
seas.list <- c('annual_quantile_975','annual_quantile_990','annual_quantile_996',
               'annual_quantile_004','annual_quantile_010','annual_quantile_025')
type <- 'annual_quantiles'

regions <- c('van_coastal_health','bella_health','northeast','willow_road')

for (region in regions) {
  for (interval in intervals) {
    for (v in seq_along(var.list)) {
      var.name <- var.list[v]
      season <- seas.list[v]
      create.ensemble.annual.climatologies(var.name,type,season,interval,gcm.list,region) 
    } 
  }
}
}

if (1==0) {
intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
var.list <- c('pr','pr','pr')
seas.list <- c('maximum','minimum','standard_deviation')
type <- 'annual'

regions <- c('van_coastal_health','bella_health','northeast','willow_road')

for (region in regions) {
  for (interval in intervals) {
    for (v in seq_along(var.list)) {
      var.name <- var.list[v]
      season <- seas.list[v]
      create.ensemble.annual.climatologies(var.name,type,season,interval,gcm.list,region) 
    } 
  }
}
}
