##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

sub.by.time <- function(input.file,interval,read.dir,write.dir) {

  nc <- nc_open(paste(read.dir,input.file,sep=''))
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
  system(paste('ncks --overwrite -d time,',st,',',en,' ',read.dir,input.file,' ',write.dir,sub.file,sep=''))
  nc_close(nc)
  rv <- paste0(write.dir,sub.file)
  return(rv)
}


run.annual.climatologies <- function(gcm,scenario) {

  intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')

  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  print(gcm)
  read.dir <- paste(proj.dir,gcm,'/',scenario,'/annual/',sep='')
  write.dir  <- paste0(read.dir,'climatologies/')
  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  ##-----------------------------------------------------------------------------
  ##Precipitation     
if (1==0) {
  all.files <- list.files(path=read.dir,pattern='pr_annual_total')
  ann.file <- all.files[grep(scenario,all.files)]

  avg.file <- gsub(pattern='pr_annual_total',replacement='pr_average_annual_total_climatology',ann.file)  
  max.file <- gsub(pattern='pr_annual_total',replacement='pr_maximum_annual_total_climatology',ann.file)  
  min.file <- gsub(pattern='pr_annual_total',replacement='pr_minimum_annual_total_climatology',ann.file)  
  sd.file  <- gsub(pattern='pr_annual_total',replacement='pr_standard_deviation_annual_total_climatology',ann.file)  
  for (interval in intervals) {        
       print(paste0('pr ',interval))
       sub.file <- sub.by.time(ann.file,interval=interval,read.dir,write.dir)
       system(paste('cdo -s -O timmean ',sub.file,' ',write.dir,gsub(pattern='1951-2100',replacement=interval,avg.file),sep=''))
       system(paste('cdo -s -O timmax ',sub.file,' ',write.dir,gsub(pattern='1951-2100',replacement=interval,max.file),sep=''))
       system(paste('cdo -s -O timmin ',sub.file,' ',write.dir,gsub(pattern='1951-2100',replacement=interval,min.file),sep=''))
       system(paste('cdo -s -O timstd ',sub.file,' ',write.dir,gsub(pattern='1951-2100',replacement=interval,sd.file),sep=''))
     system(paste('rm ',sub.file,sep=''))         
  }

  ##------------------------------------------------------------------------------
  ##Maximum Temperature
  all.files <- list.files(path=read.dir,pattern='tasmax_annual_average')
  ann.file <- all.files[grep(scenario,all.files)]

  avg.file <- gsub(pattern='tasmax_annual_average',replacement='tasmax_average_annual_climatology',ann.file)  
  for (interval in intervals) {     
       print(paste0('tx ',interval))   
       sub.file <- sub.by.time(ann.file,interval=interval,read.dir,write.dir)
       system(paste('cdo -s -O timmean ',sub.file,' ',write.dir,gsub(pattern='1951-2100',replacement=interval,avg.file),sep=''))
     system(paste('rm ',sub.file,sep=''))         
  }
  ##------------------------------------------------------------------------------
  ##Minimum Temperature
  all.files <- list.files(path=read.dir,pattern='tasmin_annual_average')
  ann.file <- all.files[grep(scenario,all.files)]

  avg.file <- gsub(pattern='tasmin_annual_average',replacement='tasmin_average_annual_climatology',ann.file)  
  for (interval in intervals) {        
       print(paste0('tn ',interval))
       sub.file <- sub.by.time(ann.file,interval=interval,read.dir,write.dir)
       system(paste('cdo -s -O timmean ',sub.file,' ',write.dir,gsub(pattern='1951-2100',replacement=interval,avg.file),sep=''))
     system(paste('rm ',sub.file,sep=''))         
  }    
}
  ##------------------------------------------------------------------------------
  ##Average Temperature
  all.tasmax.files <- list.files(path=write.dir,pattern='tasmax_average_annual_climatology')
  tasmax.files <- all.tasmax.files[grep(scenario,all.tasmax.files)]
  all.tasmin.files <- list.files(path=write.dir,pattern='tasmin_average_annual_climatology')
  tasmin.files <- all.tasmin.files[grep(scenario,all.tasmin.files)]

  for (interval in intervals) {        
       tasmax.file <- tasmax.files[grep(interval,tasmax.files)]
       tasmin.file <- tasmin.files[grep(interval,tasmin.files)]
       print(paste0('tas ',interval))
       tas.file <- gsub(pattern='tasmax',replacement='tas',tasmax.file)  
       print(tas.file)
       system(paste0('cdo -s -O add ',write.dir,tasmax.file,' ',write.dir,tasmin.file,' ',write.dir,'sub.nc'))
       system(paste0('cdo -s -O divc,2 ',write.dir,'sub.nc ',write.dir,tas.file))
       system(paste('rm ',write.dir,'sub.nc',sep=''))         
       system(paste0('ncrename -v tasmax,tas ',write.dir,tas.file))

  }    

}

run.seas.mon.temperature <- function(type,scenario='rcp85') {
  ##------------------------------------------------------------------------------
  intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
  gcms <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G','HadGEM2-CC',
            'HadGEM2-ES','inmcm4','MRI-CGCM3','MIROC5','MPI-ESM-LR','CCSM4')

  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  for (gcm in gcms) {
    print(gcm)
    read.dir <- paste(proj.dir,gcm,'/',scenario,'/',type,'/',sep='')
    write.dir  <- paste0(read.dir,'climatologies/')
    if (!file.exists(write.dir)) {
        dir.create(write.dir,recursive=T)
    }
  
    ##Average Temperature
    all.tasmax.files <- list.files(path=write.dir,pattern=paste0('tasmax_',type,'_average_climatology'))
    tasmax.files <- all.tasmax.files[grep(scenario,all.tasmax.files)]
    all.tasmin.files <- list.files(path=write.dir,pattern=paste0('tasmin_',type,'_average_climatology'))
    tasmin.files <- all.tasmin.files[grep(scenario,all.tasmin.files)]

    for (interval in intervals) {        
       tasmax.file <- tasmax.files[grep(interval,tasmax.files)]
       tasmin.file <- tasmin.files[grep(interval,tasmin.files)]
       print(paste0('tas ',interval))
       tas.file <- gsub(pattern='tasmax',replacement='tas',tasmax.file)  
       print(tas.file)

       system(paste0('cdo -s -O add ',write.dir,tasmax.file,' ',write.dir,tasmin.file,' ',write.dir,'sub.nc'))
       system(paste0('cdo -s -O divc,2 ',write.dir,'sub.nc ',write.dir,tas.file))
       system(paste('rm ',write.dir,'sub.nc',sep=''))         
       system(paste0('ncrename -v tasmax,tas ',write.dir,tas.file))
    }
  }    
}

run.climdex.climatologies <- function(var.list,gcm,scenario='rcp85') {

  intervals <- c('1971-2000','2011-2040','2041-2070','2071-2099')

  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  print(gcm)
  read.dir <- paste(proj.dir,gcm,'/',scenario,'/climdex/',sep='')
  write.dir  <- paste0(read.dir,'climatologies/')
  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  for (var.name in var.list) {
      print(var.name)
      seas.fx <- get.seas.fx(var.name)
      clim.file <- list.files(path=read.dir,pattern=paste0('^',var.name))
      print(clim.file)
      monthly <- grepl('_mon_',clim.file)
      if (monthly) {      
         pattern <- paste0(var.name,'_mon')
         mon.replacement <- paste(var.name,'_monthly_climatology',sep='')  
         month.new <- gsub(pattern=pattern,replacement=mon.replacement,clim.file)
         ann.replacement <- paste(var.name,'_annual_climatology',sep='')  
         annual.new <- gsub(pattern=pattern,replacement=ann.replacement,clim.file)
         seas.replacement <- paste(var.name,'_seasonal_climatology',sep='')  
         seasonal.new <- gsub(pattern=pattern,replacement=seas.replacement,clim.file)
              
         for (interval in intervals) {     
           write.month <- gsub(pattern='1951-2099',replacement=interval,month.new)
           sub.file <- sub.by.time(clim.file,interval=interval,read.dir,write.dir)
           system(paste('cdo -s -O ymonmean ',sub.file,' ',write.dir,write.month,sep=''))
           system(paste('rm ',sub.file,sep=''))         

           write.seas <- gsub(pattern='1951-2099',replacement=interval,seasonal.new)
           sub.file <- sub.by.time(clim.file,interval=interval,read.dir,write.dir)
           system(paste('cdo -s -O seas',seas.fx,' ',sub.file,' ',write.dir,'tmp.nc',sep=''))
           system(paste('cdo -s -O yseasmean ',write.dir,'tmp.nc ',write.dir,write.seas,sep=''))
           system(paste('rm ',sub.file,sep=''))         

           write.year <- gsub(pattern='1951-2099',replacement=interval,annual.new)
           sub.file <- sub.by.time(clim.file,interval=interval,read.dir,write.dir)
           system(paste('cdo -s -O year',seas.fx,' ',sub.file,' ',write.dir,'tmp.nc',sep=''))
           system(paste('cdo -s -O timmean ',write.dir,'tmp.nc ',write.dir,write.year,sep=''))
           system(paste('rm ',sub.file,sep=''))         
           system(paste('rm ',write.dir,'tmp.nc',sep=''))         
         }
      } else {
         pattern <- paste0(var.name,'_ann')
         replacement <- paste(var.name,'_annual_climatology',sep='')  
         annual.new <- gsub(pattern=pattern,replacement=replacement,clim.file)
             
         for (interval in intervals) {     
           write.year <- gsub(pattern='1951-2099',replacement=interval,annual.new)
           sub.file <- sub.by.time(clim.file,interval=interval,read.dir,write.dir)
           system(paste('cdo -s -O timmean ',sub.file,' ',write.dir,write.year,sep=''))
           system(paste('rm ',sub.file,sep=''))         
         }
      }
  } 
}

run.drought.climatologies <- function(gcm,scenario) {

  intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')

  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  print(gcm)
  read.dir <- paste(proj.dir,gcm,'/',scenario,'/climdex/',sep='')
  write.dir  <- paste0(read.dir,'climatologies/')
  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  ##-----------------------------------------------------------------------------
  ##Precipitation     
  clim.files <- list.files(path=read.dir,pattern='cddETCCDI')
  clim.file <- clim.files[grep(scenario,clim.files)]

  cdd.30.file <- gsub(pattern='cddETCCDI',replacement='cdd30ETCCDI_annual_climatology',clim.file)  
  cdd.90.file <- gsub(pattern='cddETCCDI',replacement='cdd90ETCCDI_annual_climatology',clim.file)  
  cdd.max.file <- gsub(pattern='cddETCCDI',replacement='cddmaxETCCDI_annual_climatology',clim.file)  

  for (interval in intervals) {        
     print(paste0('cdd ',interval))
     sub.file <- sub.by.time(clim.file,interval=interval,read.dir,write.dir)
     system(paste('cdo -s -O timsum -gec,30 ',sub.file,' ',write.dir,gsub(pattern='1951-2100',replacement=interval,cdd.30.file),sep=''))
##     system(paste('cdo -s -O timmax ',sub.file,' ',write.dir,gsub(pattern='1951-2099',replacement=interval,cdd.max.file),sep=''))
##     system(paste('cdo -s -O timpctl,90 ',sub.file,' -timmin ',sub.file,' -timmax ',sub.file,' ',
##                   write.dir,gsub(pattern='1951-2099',replacement=interval,cdd.90.file),sep=''))
     system(paste('rm ',sub.file,sep=''))         
  }
}
  
run.standard.climatologies <- function(var.list,gcm,scenario='rcp85',type,grep.fx,clim.fx) {

  intervals <- c('1971-2000','2011-2040','2041-2070','2071-2099')

  proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

  print(gcm)
  read.dir <- paste(proj.dir,gcm,'/',scenario,'/',type,'/',sep='')
  write.dir  <- paste0(read.dir,'climatologies/')
  if (!file.exists(write.dir)) {
     dir.create(write.dir,recursive=T)
  }

  for (var.name in var.list) {
      print(var.name)
      if (type=='annual_quantiles') {
         grep.name <- grep.fx 
      } else {
        grep.name <- grep.fx(var.name)
      }
      all.files <- list.files(path=read.dir,pattern=paste0(var.name,'_',grep.name))
      mon.file <- all.files[grep(scenario,all.files)]

      pattern <- paste0(var.name,'_',grep.name)
      if (type=='seasonal' | type=='monthly') {
         grep.name <- paste0(grep.name,'_average')
      }
      replacement <- paste(var.name,'_',grep.name,'_climatology',sep='')  
      file.new <- gsub(pattern=pattern,replacement=replacement,mon.file)
     
    for (interval in intervals) {     
      write.new <- gsub(pattern='1951-2099',replacement=interval,file.new)
      sub.file <- sub.by.time(mon.file,interval=interval,read.dir,write.dir)
      system(paste('cdo -s -O ',clim.fx,' ',sub.file,' ',write.dir,write.new,sep=''))
      system(paste('rm ',sub.file,sep=''))         

    }
  } 
}

get.seas.fx <- function(var.name) {
  fx <- switch(var.name,
               rx1dayETCCDI='max',
               rx2dayETCCDI='max',
               rx5dayETCCDI='max',
               txxETCCDI='max',
               tnxETCCDI='max',
               txnETCCDI='min',
               tnnETCCDI='min',
               tn10pETCCDI='mean',
               tx10pETCCDI='mean',
               tn90pETCCDI='mean',
               tx90pETCCDI='mean',
               dtrETCCDI='mean')
  if (is.null(fx))
    fx <- mean
  return(fx)
}

grep.ext <- function(var.name) {
              rv <- switch(var.name,
                          pr='annual_maximum',
                          tasmax='annual_maximum',
                          tasmin='annual_minimum') 
} 
grep.seas <- function(var.name) {return('seasonal')}
grep.mon <- function(var.name) {return('monthly')}
grep.dd <-  function(var.name) {return('annual')}

##************************************************************************

gcm <- 'MRI-CGCM3'
scenario <- 'rcp85'
type <- 'monthly'

##  args <- commandArgs(trailingOnly=TRUE)
##  for(i in 1:length(args)){
##      eval(parse(text=args[[i]]))
##  }

run.seas.mon.temperature(type,scenario)
browser()
##Annual
if (type=='annual') {
  run.annual.climatologies(gcm,scenario)
}
##-----------------------------------------
##Extreme Values
if (type=='annual_extremes') {
  var.list <- c('pr','tasmax','tasmin')
  run.standard.climatologies(var.list,gcm,scenario='rcp85',type='annual_extremes',
                             grep.fx=grep.ext,clim.fx='timmean')
}

##Monthly
if (type=='monthly') {
  var.list <- c('pr','tasmax','tasmin')
  run.standard.climatologies(var.list,gcm,scenario='rcp85',type='monthly',
                             grep.fx=grep.mon,clim.fx='ymonmean')
}
##Seasonal
if (type=='seasonal') {
  var.list <- c('pr','tasmax','tasmin')
  run.standard.climatologies(var.list,gcm,scenario='rcp85',type='seasonal',
                             grep.fx=grep.seas,clim.fx='yseasmean')
}

##------------------------------------------
##Degree Days
if (type=='degree_days') {
  var.list <- c('cdd','fdd','gdd','hdd')
  run.standard.climatologies(var.list,gcm,scenario='rcp85',type='degree_days',
                             grep.fx=grep.dd,clim.fx='timmean')
}
##------------------------------------------
##Climdex              
if (type=='climdex') {
###'sdiiETCCDI','trETCCDI',

##  climdex.names <- c('gslETCCDI','fdETCCDI','suETCCDI','su30ETCCDI','idETCCDI',
##                     'rx1dayETCCDI','rx2dayETCCDI', 'rx5dayETCCDI',
##                     'txxETCCDI','tnxETCCDI','txnETCCDI','tnnETCCDI',                    
##                     'r10mmETCCDI','r20mmETCCDI','cddETCCDI','cwdETCCDI',
##                     'r95pETCCDI','r99pETCCDI','prcptotETCCDI','r95daysETCCDI','r99daysETCCDI',
##                     'dtrETCCDI')
climdex.names <- c('r95daysETCCDI')

##gcms <- c('ACCESS1-0','CanESM2','CCSM4','CSIRO-Mk3-6-0','CNRM-CM5','inmcm4')
##for (gcm in gcms) {
    run.climdex.climatologies(climdex.names,gcm,scenario='rcp85')
##}

}

##------------------------------------------
##CDD Drought
if (type=='drought') {
  gcms <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G','HadGEM2-CC',
            'HadGEM2-ES','inmcm4','MRI-CGCM3','MIROC5','MPI-ESM-LR','CCSM4')
 ## gcms <- c('HadGEM2-CC','HadGEM2-ES')
  for (gcm in gcms) {
      run.drought.climatologies(gcm,scenario='rcp85')
  }
}

##------------------------------------------
##Annual Quantiles
if (type=='annual_quantiles') {
  gcms <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G','inmcm4','MRI-CGCM3',
            'MIROC5','MPI-ESM-LR','CCSM4')
  gcms <- c('HadGEM2-CC','HadGEM2-ES')
  for (gcm in gcms) {
    run.standard.climatologies('tasmax',gcm,scenario='rcp85',type='annual_quantiles',grep.fx='annual_quantile_975',clim.fx='timmean')
    run.standard.climatologies('tasmax',gcm,scenario='rcp85',type='annual_quantiles',grep.fx='annual_quantile_990',clim.fx='timmean')
    run.standard.climatologies('tasmax',gcm,scenario='rcp85',type='annual_quantiles',grep.fx='annual_quantile_996',clim.fx='timmean')
    run.standard.climatologies('tasmin',gcm,scenario='rcp85',type='annual_quantiles',grep.fx='annual_quantile_004',clim.fx='timmean')
    run.standard.climatologies('tasmin',gcm,scenario='rcp85',type='annual_quantiles',grep.fx='annual_quantile_010',clim.fx='timmean')
    run.standard.climatologies('tasmin',gcm,scenario='rcp85',type='annual_quantiles',grep.fx='annual_quantile_025',clim.fx='timmean')
  }
}
