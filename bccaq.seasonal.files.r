##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

##RCP26
## 
## 
## 
##
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


avg.sub.time <- function(day.file,seas.file,var.name,interval,read.dir,write.dir) {

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
  
  st <- head(grep(yrs[1],years),1)-1
  en <- tail(grep(yrs[2],years),1)-1
  sub.file <- 'sub.nc'
  system(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,day.file,' ',write.dir,sub.file,sep=''))
  print(seas.file)   
  ##For precip
  if (var.name=='pr') {
    ##Seasonal/Monthly
    system(paste('cdo -s -O seassum ',write.dir,sub.file,' ',write.dir,seas.file,sep=''))
    ##system(paste('cdo -s -O seassum ',write.dir,sub.file,' ',write.dir,'tmp.nc',sep=''))
    ##system(paste('cdo -s -O yseasmean ',write.dir,'tmp.nc ',write.dir,seas.file,sep=''))
    ##Annual
    ##system(paste('cdo -s -O yearsum ',write.dir,sub.file,' ',write.dir,seas.file,sep=''))
    ##system(paste('cdo -s -O yearsum ',write.dir,sub.file,' ',write.dir,'tmp.nc',sep=''))
    ##system(paste('cdo -s -O timmean ',write.dir,'tmp.nc ',write.dir,seas.file,sep=''))

    ##system(paste('rm ',write.dir,'tmp.nc',sep=''))
    system(paste('rm ',write.dir,'sub.nc',sep=''))
  } else {
    ##For temperature  
    ##Seasonal/Monthly
    system(paste('cdo -s -O seasmean ',write.dir,sub.file,' ',write.dir,seas.file,sep=''))
    ##system(paste('cdo -s -O yseasmean ',write.dir,sub.file,' ',write.dir,seas.file,sep=''))
    ##Annual
    ##system(paste('cdo -s -O yearmean ',write.dir,sub.file,' ',write.dir,seas.file,sep=''))
    ##system(paste('cdo -s -O yearmean ',write.dir,sub.file,' ',write.dir,'tmp.nc',sep=''))
    ##system(paste('cdo -s -O timmean ',write.dir,'tmp.nc ',write.dir,seas.file,sep=''))
    ##system(paste('rm ',write.dir,'tmp.nc',sep=''))
    system(paste('rm ',write.dir,'sub.nc',sep=''))
  }  
  nc_close(nc)
}

run.metro.van.seasons <- function() {
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

  ds.type <- 'bccaq'
  var.list <- c('pr','tasmax','tasmin') ##c('pr','tasmax','tasmin') ##'s30' ##
  past.int <- '1951-2000'
  ##proj.int <- '2011-2040'
  scenario <- 'rcp85'
  
  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  for (gcm in gcm.list) {
    print(gcm)
    read.dir <- paste(proj.dir,gcm,'/',sep='')
    write.dir  <- paste(proj.dir,scenario,'/seasonal/',gcm,'/',sep='')
    
    if (!file.exists(write.dir))
      dir.create(write.dir)
    for (var.name in var.list) {
      print(var.name)
      all.files <- list.files(path=read.dir,pattern=paste(var.name,'_day',sep=''))
      ##all.files <- list.files(path=read.dir,pattern=paste(var.name,'_gcm_prism',sep=''))
      scenario.files <- all.files[grep(scenario,all.files)]
      past.file <- scenario.files[grep('1951-2000',scenario.files)]
      proj.file <- scenario.files[grep('2001-2100',scenario.files)]
      ##past.file <- proj.file <- scenario.files
##browser()
      pattern <- paste(var.name,'_day',sep='')
      ##pattern <- paste(var.name,'_gcm',sep='')
      replacement <- paste(var.name,'_seas',sep='')  
      var.past <- gsub(pattern=pattern,replacement=replacement,past.file)
      var.proj <- gsub(pattern=pattern,replacement=replacement,proj.file)

      ##write.past <- gsub(pattern='19500101-21001231',replacement=past.int,var.past)
      ##write.proj <- gsub(pattern='19500101-21001231',replacement=proj.int,var.proj)
     
      write.past <- gsub(pattern='1951-2000',replacement=past.int,var.past)
      avg.sub.time(past.file,write.past,var.name,interval=past.int,read.dir,write.dir)

      write.proj <- gsub(pattern='2001-2100',replacement='2001-2100',var.proj)
      avg.sub.time(proj.file,write.proj,var.name,interval='2001-2100',read.dir,write.dir)

##      write.proj <- gsub(pattern='2001-2100',replacement='2011-2040',var.proj)
##      avg.sub.time(proj.file,write.proj,var.name,interval='2011-2040',read.dir,write.dir)

##      write.proj <- gsub(pattern='2001-2100',replacement='2041-2070',var.proj)
##      avg.sub.time(proj.file,write.proj,var.name,interval='2041-2070',read.dir,write.dir)

##      write.proj <- gsub(pattern='2001-2100',replacement='2071-2100',var.proj)
##      avg.sub.time(proj.file,write.proj,var.name,interval='2071-2100',read.dir,write.dir)

    } 

  }
}


make.tas.files <- function() {

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
  
  ds.type <- 'bccaq'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  
  proj.dir <- '/storage/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_bc_subset/rcp85/' ##'/home/data/scratch/ssobie/bccaq_gcm_okanagan_subset/'
  
  for (gcm in gcm.list) {
    print(gcm)
    read.dir <- paste(proj.dir,gcm,'/',sep='')
    write.dir  <- read.dir
    
    tx.files <- list.files(path=read.dir,pattern='tasmax_day')
    tx.past.file <- tx.files[grep(past.int,tx.files)]
    tx.proj.file <- tx.files[grep(proj.int,tx.files)]

    tn.files <- list.files(path=read.dir,pattern='tasmin_day')
    tn.past.file <- tn.files[grep(past.int,tn.files)]
    tn.proj.file <- tn.files[grep(proj.int,tn.files)]
    
    tas.past <- gsub(pattern='tasmax',replacement='tas',tx.past.file)
    tas.proj <- gsub(pattern='tasmax',replacement='tas',tx.proj.file)

    system(paste('cdo -s add ',read.dir,tx.past.file,' ',read.dir,tn.past.file,' ',write.dir,'tmp.nc',sep=''))
    system(paste('cdo -s divc,2 ',write.dir,'tmp.nc ',write.dir,tas.past,sep=''))
    system(paste('ncrename -v tasmax,tas ',write.dir,tas.past,sep=''))
    
    system(paste('cdo -s add ',read.dir,tx.proj.file,' ',read.dir,tn.proj.file,' ',write.dir,'tmp.nc',sep=''))
    system(paste('cdo -s divc,2 ',write.dir,'tmp.nc ',write.dir,tas.proj,sep=''))
    system(paste('ncrename -v tasmax,tas ',write.dir,tas.proj,sep=''))
    
    system(paste('rm ',write.dir,'tmp.nc ',sep=''))
    ##system(paste('rm ',write.dir,'*anoms_interp*',sep=''))


    
  }
}



rename.tas.files <- function() {

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

  ds.type <- 'bccaq'
  past.int <- '1951-2000'
  proj.int <- '2001-2100'
  
  proj.dir <-  '/home/data/scratch/ssobie/bccaq_gcm_okanagan_subset/'
  
  for (gcm in gcm.list) {
    print(gcm)
    read.dir <- paste(proj.dir,gcm,'/',sep='')
    write.dir  <- read.dir
    
    tas.files <- list.files(path=read.dir,pattern='tas_gcm_prism')
    for (file in tas.files) {
      print(file)
      system(paste('ncrename -v tasmax,tas ',read.dir,file,sep=''))
    }

  }
}


run.metro.van.seasons()