##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

ann.sub.time <- function(ann.file,clim.file,var.name,interval,read.dir,write.dir) {

  nc <- nc_open(paste(read.dir,ann.file,sep=''))
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
  system(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,ann.file,' ',write.dir,sub.file,sep=''))
  ##For precip
    ##Annual
    system(paste('cdo -s -O timmean ',write.dir,sub.file,' ',write.dir,clim.file,sep=''))
    system(paste('rm ',write.dir,'sub.nc',sep=''))
  nc_close(nc)
}

mon.sub.time <- function(mon.file,clim.file,var.name,interval,read.dir,write.dir) {

  nc <- nc_open(paste(read.dir,mon.file,sep=''))
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
  system(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,mon.file,' ',write.dir,sub.file,sep=''))
  ##For precip
    ##Monthly
    system(paste('cdo -s -O ymonmean ',write.dir,sub.file,' ',write.dir,clim.file,sep=''))
    system(paste('rm ',write.dir,'sub.nc',sep=''))  
  nc_close(nc)
}


run.annual.climatologies <- function() {

  var.list <- c('pr','tasmax','tasmin')
  past.int <- '1971-2000'
  scenario <- 'rcp85'
  
  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  gcm.list <- 'CanESM2'

  for (gcm in gcm.list) {
    print(gcm)
    read.dir <- paste(proj.dir,gcm,'/',scenario,'/annual/',sep='')
    write.dir  <- read.dir
    
    for (var.name in var.list) {
      print(var.name)
      grep.name <- switch(var.name,
                          pr='_annual_total_',
                          tasmax='_annual_average_',
                          tasmin='_annual_average_')

      all.files <- list.files(path=read.dir,pattern=paste(var.name,grep.name,sep=''))
      ann.file <- all.files[grep(scenario,all.files)]

      pattern <- paste(var.name,grep.name,sep='')
      replacement <- paste(var.name,'_annual_climatology_',sep='')  
      file.new <- gsub(pattern=pattern,replacement=replacement,ann.file)
     
      write.past <- gsub(pattern='1951-2100',replacement=past.int,file.new)
      ann.sub.time(ann.file,write.past,var.name,interval=past.int,read.dir,write.dir)

      write.proj <- gsub(pattern='1951-2100',replacement='2011-2040',file.new)
      ann.sub.time(ann.file,write.proj,var.name,interval='2011-2040',read.dir,write.dir)

      write.proj <- gsub(pattern='1951-2100',replacement='2041-2070',file.new)
      ann.sub.time(ann.file,write.proj,var.name,interval='2041-2070',read.dir,write.dir)

      write.proj <- gsub(pattern='1951-2100',replacement='2071-2100',file.new)
      ann.sub.time(ann.file,write.proj,var.name,interval='2071-2100',read.dir,write.dir)

    } 

  }
}


run.monthly.climatologies <- function() {

  var.list <- c('pr','tasmax','tasmin')
  past.int <- '1971-2000'
  scenario <- 'rcp85'
  grep.name <- 'monthly'  
  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  gcm.list <- 'CanESM2'

  for (gcm in gcm.list) {
    print(gcm)
    read.dir <- paste(proj.dir,gcm,'/',scenario,'/monthly/',sep='')
    write.dir  <- read.dir
    
    for (var.name in var.list) {
      print(var.name)

      all.files <- list.files(path=read.dir,pattern=paste(var.name,'_monthly_',sep=''))
      mon.file <- all.files[grep(scenario,all.files)]

      pattern <- paste(var.name,'_monthly_',sep='')
      replacement <- paste(var.name,'_monthly_climatology_',sep='')  
      file.new <- gsub(pattern=pattern,replacement=replacement,mon.file)
     
      write.past <- gsub(pattern='1951-2100',replacement=past.int,file.new)
      mon.sub.time(mon.file,write.past,var.name,interval=past.int,read.dir,write.dir)

      write.proj <- gsub(pattern='1951-2100',replacement='2011-2040',file.new)
      mon.sub.time(mon.file,write.proj,var.name,interval='2011-2040',read.dir,write.dir)

      write.proj <- gsub(pattern='1951-2100',replacement='2041-2070',file.new)
      mon.sub.time(mon.file,write.proj,var.name,interval='2041-2070',read.dir,write.dir)

      write.proj <- gsub(pattern='1951-2100',replacement='2071-2100',file.new)
      mon.sub.time(mon.file,write.proj,var.name,interval='2071-2100',read.dir,write.dir)

    } 

  }
}


