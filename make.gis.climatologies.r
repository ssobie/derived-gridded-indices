##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)
library(raster)
options(warn=2)
##-----------------------------------------------------------------

make.netcdf.file <- function(file.name,var.name,clim.data) {

    nc <- nc_open(file.name,write=TRUE)
    if (grepl('(cdd90|cddmax)',var.name)) {
        var.name <- 'cddETCCDI'
    }
    ncvar_put(nc,var.name,clim.data)  
    ncatt_put(nc,varid=0,attname='history',attval='')
    ncatt_put(nc,varid=0,attname='driving_experiment',
                      attval='PCIC 12 Ensemble Average, historical+rcp85')
    ncatt_put(nc,varid=0,attname='driving_model_id',
                      attval='PCIC 12 Ensemble Average')
    nc_close(nc)    
}

##--------------------------------------------------------------

##Annual Climatologies
create.ensemble.annual.climatologies <- function(var.name,type,season,interval,gcm.list,region) {

  scenario <- 'rcp85'
 
  if (region == 'bc') {
    proj.dir <-   '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/' ##For BC
  } else {
    proj.dir <-   paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/')
    ##proj.dir <-   paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/')
  }

  write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/',scenario,'/',type,'/ENSEMBLE/')

  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }
  gx <- length(gcm.list)

  if (region=='bc') {
    copy.dir <- paste0(proj.dir,'ACCESS1-0/rcp85/',type,'/climatologies/') ### For BC
    if (type=='return_periods') {
           copy.dir <- paste0(proj.dir,'ACCESS1-0/rcp85/',type,'/ACCESS1-0/') ### For BC Return Periods
    }
  } else {
    copy.dir <- paste0(proj.dir,'rcp85/',type,'/ACCESS1-0/')
  }
  ##New files     
  copy.files <- list.files(path=copy.dir,pattern=interval)

  var.files <- copy.files[grep(paste0('^',var.name,'_'),copy.files)]
  copy.file <- var.files[grep(paste0(season,'_'),var.files)]

  new.file <- gsub('ACCESS1-0','ENSEMBLE',copy.file)
  if (region=='bc') {
     if (type=='return_periods') {
        new.file <- gsub(var.name,paste0(var.name,'_bc'),new.file)
        new.file <- gsub('BCCAQ2','climatology_BCCAQ2',new.file)
     }
  }

  file.copy(from=paste0(copy.dir,copy.file),
            to=paste0(write.dir,new.file),overwrite=T)
              
  file.split <- strsplit(copy.file,'_')[[1]]

  data.ens <- c()

  for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    if (region=='bc') {
       read.dir <- paste(proj.dir,'/',gcm,'/rcp85/',type,'/climatologies/',sep='') ##For BC
       if (type=='return_periods') {
          read.dir <- paste0(proj.dir,gcm,'/rcp85/',type,'/',gcm,'/') ### For BC Return Periods
       }
    } else {
       read.dir <- paste(proj.dir,'/rcp85/',type,'/',gcm,'/',sep='')
    }
    all.files <- list.files(path=read.dir,pattern=interval)
    var.files <- all.files[grep(paste0('^',var.name,'_'),all.files)]
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
  ##Return period names
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

##------------------------------------------------------------------------

##Seasonal and Monthly climatologies
create.ensemble.seas.mon.climatologies <- function(var.name,type,season,interval,gcm.list,region) {

  scenario <- 'rcp85'
  if (region == 'north_america') {
    proj.dir <-   '/storage/data/climate/downscale/CMIP5/building_code/' ##For BC
  } else if (region == 'bc') {
    proj.dir <-   '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/' ##For BC
  } else {
    proj.dir <-   paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/')
  }

  if (region =='north_america') {
    write.dir <- '/storage/data/climate/downscale/CMIP5/building_code/ENSEMBLE/'
  } else {
    write.dir  <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/',region,'/',scenario,'/',type,'/ENSEMBLE/')
  }
  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  seas.ix <- switch(season,
                    winter=1,spring=2,summer=3,fall=4,
                    january=1,february=2,march=3,april=4,may=5,june=6,july=7,august=8,september=9,october=10,november=11,december=12)
  gx <- length(gcm.list)

  if (region == 'north_america') {
     copy.dir <- paste0(proj.dir,'ACCESS1-0/climatologies/') ##For BC
  } else if (region == 'bc') {
     copy.dir <- paste0(proj.dir,'ACCESS1-0/rcp85/',type,'/climatologies/') ##For BC
  } else {
     copy.dir <- paste0(proj.dir,'rcp85/',type,'/ACCESS1-0/')
  }
  ##New files 
  copy.files <- list.files(path=copy.dir,pattern=interval)
  var.files <- copy.files[grep(paste0(var.name,'_'),copy.files)]

  if (grepl('(seas|mon)',type)) {
    copy.file <- var.files[grep(type,var.files)]
    inter.file <- gsub('ACCESS1-0','ENSEMBLE',copy.file)
    new.file <- gsub(type,season,inter.file)      
    if (region=='north_america') {                  
      new.file <- gsub('monthly',season,inter.file)
    }
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
    ##print(gcm)
    if (region=='north_america') {
       read.dir <- paste(proj.dir,'/',gcm,'/climatologies/',sep='') ##For BC
    } else if (region=='bc') {
       read.dir <- paste(proj.dir,'/',gcm,'/rcp85/',type,'/climatologies/',sep='') ##For BC
    } else {
       read.dir <- paste(proj.dir,'/rcp85/',type,'/',gcm,'/',sep='')
    }
    all.files <- list.files(path=read.dir,pattern=interval)
    var.files <- all.files[grep(paste0(var.name,'_'),all.files)]
    if (grepl('(seas|mon)',type)) {
      seas.file <- var.files[grep(type,var.files)]
    } else { ##For climdex files
      seas.file  <- var.files[grep('seasonal',var.files)]
    }
    print(seas.file)
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

##--------------------------------------------------------------

##Compute anomalies (absolute or percent) from the climatologies 
create.anomalies.from.climatologies <- function(var.name,type,season,region,proj.int,
                                                past.int='1971-2000') {
  base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/assessment_subsets/'
  read.dir  <- paste0(base.dir,region,'/rcp85/',type,'/ENSEMBLE/')
  write.dir  <- paste0(base.dir,region,'/rcp85/',type,'/ENSEMBLE/anomalies/')
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


##*************************************************************************
##*************************************************************************
gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0',
              'GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4',
              'MIROC5','MPI-ESM-LR','MRI-CGCM3')

regions <- 'bc' ##c('central','kootenays')
intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')

for (region in regions) {
  for (interval in intervals) {
if (1==1) {
    ##Degree Day Climatologies
    var.list <- c('cdd','fdd','gdd','hdd')
    for (var.name in var.list) {
##        create.ensemble.annual.climatologies(var.name,type='degree_days',season='annual',interval,gcm.list,region) 
    }

    var.list <- 'wetbulb' ##'tas' ##c('pr','tasmax','tasmin','tas')

    ##Annual Climatologies
    for (var.name in var.list) {
##        create.ensemble.annual.climatologies(var.name,type='annual',season='annual',interval,gcm.list,region) 
    } 

    ##Monthly Climatologies
    seasons <- c('january','february','march','april','may','june','july','august','september','october','november','december')
    for (var.name in var.list) {
      for (season in seasons) {
##        create.ensemble.seas.mon.climatologies(var.name,type='monthly_001',season,interval,gcm.list,region) 
      }
    }

    ##Seasonal Climatologies
    seasons <- c('winter','spring','summer','fall')
    for (var.name in var.list) {
      for (season in seasons) {
##        create.ensemble.seas.mon.climatologies(var.name,type='seasonal',season,interval,gcm.list,region) 
      }
    }

    ##Return Periods
    var.list <- c('pr','tasmax','tasmin')
    for (var.name in var.list) {
        create.ensemble.annual.climatologies(var.name,type='return_periods',season='RP20',interval,gcm.list,region) 

    } 
} 
if (1==0) {
    ##----------------------------------------------------
    ##Annual Climdex
    var.list <-  c('fdETCCDI','suETCCDI','su30ETCCDI','idETCCDI','trETCCDI','gslETCCDI',
                   'txxETCCDI','txnETCCDI','tnnETCCDI','tnxETCCDI','dtrETCCDI',
                   'rx1dayETCCDI','rx5dayETCCDI',
                   'sdiiETCCDI','r10mmETCCDI','r20mmETCCDI','cwdETCCDI',
                   'cwdETCCDI','cddETCCDI','prcptotETCCDI','cdd90ETCCDI','cddmaxETCCDI',
                   'r95pETCCDI','r99pETCCDI','r95daysETCCDI','r99daysETCCDI')  
  var.list <- c('txxETCCDI','tnnETCCDI')    
    for (var.name in var.list) {
         create.ensemble.annual.climatologies(var.name,type='climdex',season='annual',interval,gcm.list,region)
    }

    ##-----------------------------------------------------
    ##Seasonal Climdex
    var.list <-  c('txxETCCDI','tnnETCCDI') ##,'txnETCCDI','tnxETCCDI','dtrETCCDI',
                   ##'rx1dayETCCDI','rx5dayETCCDI')
    seasons <- c('winter','spring','summer','fall')
    for (var.name in var.list) {
      print(var.name)
      for (season in seasons) {
        create.ensemble.seas.mon.climatologies(var.name,type='climdex',season,interval,gcm.list,region) 
      }
    }
}
  }##Intervals
}##Regions



##Annual Quantiles
if (1==0) {
intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
var.list <- c('tasmax','tasmax','tasmax',
              'tasmin','tasmin','tasmin')
seas.list <- c('annual_quantile_975','annual_quantile_990','annual_quantile_996',
               'annual_quantile_004','annual_quantile_010','annual_quantile_025')
type <- 'annual_quantiles'

##regions <- c('van_coastal_health','bella_health','northeast','willow_road')
##regions <- c('interior_health','toquaht')
##region <- 'fraser_health'
region <- 'bc'

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

regions <- 'interior_health' ##c('van_coastal_health','bella_health','northeast','willow_road')

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
