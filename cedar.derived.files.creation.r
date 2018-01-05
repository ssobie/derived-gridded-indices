##Script to calculate the degree-day indices from the BCCAQ data extracted from the large files
#and write the output to a new netcdf file

##Updated version from compute.climdex.bccaq.r
##This computes all the climdex variables

library(ncdf4)
library(PCICt)

time.component <- function(time.file,freq) {

  nc <- nc_open(time.file,write=FALSE)

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                          cal=time.calendar)
 
  past.values <- ncvar_get(nc,'time')
  full.series <- format(past.origin + past.values*86400,'%Y-%m-%d')
  full.years <- format(past.origin + past.values*86400,'%Y')
  print(range(full.years))
  years.ix <- grep('*-01-01',full.series)
  years <- past.values[years.ix]
  seas.ix <- sort(c(grep('*-01-01',full.series),
		    grep('*-04-01',full.series),
		    grep('*-08-01',full.series),
		    grep('*-12-01',full.series)))
  seasons <- past.values[seas.ix]		    
  months.ix <- grep('[0-9]{4}-[0-9]{2}-01',full.series)
  months <- past.values[months.ix]

  dates <- switch(freq,
                  Ann=years,
                  Mon=months,
		  Seas=seasons,
                  AnnClim=1,
                  SeasClim=1:4,
                  MonClim=1:12,
		  RP=1:3)

  new.int <- paste(range(full.years),collapse='-')
  print(new.int)
  dates <- as.numeric(dates)
  nc_close(nc)

  rv <- list(atts=time.atts,
             dates=dates,
	     interval=new.int)
  return(rv)
      
}

space.component <- function(space.file) {

  print(space.file)
  nc <- nc_open(space.file,write=FALSE)
  ##Attributes to retain
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')  

  lon.atts <- ncatt_get(nc,'lon')
  lat.atts <- ncatt_get(nc,'lat')
  global.atts <- ncatt_get(nc,varid=0)
  nc_close(nc)
  rv <- list(global.atts=global.atts,
             lon.atts=lon.atts,
             lat.atts=lat.atts,
             lon=lon,
             lat=lat)  
  return(rv)   
}

get.climdex.info <- function(climdex.name) {

  climdex.names <- list(climdex.fd=c('fdETCCDI','Ann','days'),
                        climdex.su=c('suETCCDI','Ann','days'),
                        climdex.id=c('idETCCDI','Ann','days'),
                        climdex.tr=c('trETCCDI','Ann','days'),
                        climdex.gsl=c('gslETCCDI','Ann','days'),
                        climdex.txx=c('txxETCCDI','Mon','degC'),
                        climdex.tnx=c('tnxETCCDI','Mon','degC'),
                        climdex.txn=c('txnETCCDI','Mon','degC'),
                        climdex.tnn=c('tnnETCCDI','Mon','degC'),
                        climdex.tn10p=c('tn10pETCCDI','Mon','days'),
                        climdex.tx10p=c('tx10pETCCDI','Mon','days'),
                        climdex.tn90p=c('tn90pETCCDI','Mon','days'),
                        climdex.tx90p=c('tx90pETCCDI','Mon','days'),
                        climdex.wsdi=c('wsdiETCCDI','Ann','days'),
                        climdex.csdi=c('csdiETCCDI','Ann','days'),
                        climdex.dtr=c('dtrETCCDI','Mon','degC'),
                        climdex.rx1day=c('rx1dayETCCDI','Mon','mm'),
                        climdex.rx2day=c('rx2dayETCCDI','Mon','mm'),
                        climdex.rx5day=c('rx5dayETCCDI','Mon','mm'),
                        climdex.sdii=c('sdiiETCCDI','Ann','mm d-1'),
                        climdex.r10mm=c('r10mmETCCDI','Ann','days'),
                        climdex.r20mm=c('r20mmETCCDI','Ann','days'),
                        climdex.cdd=c('cddETCCDI','Ann','days'),
                        climdex.cwd=c('cwdETCCDI','Ann','days'),
                        climdex.r95ptot=c('r95pETCCDI','Ann','mm'),
                        climdex.r99ptot=c('r99pETCCDI','Ann','mm'),
                        climdex.prcptot=c('prcptotETCCDI','Ann','mm'),
                        climdex.r95days=c('r95daysETCCDI','Ann','days'),
                        climdex.r99days=c('r99daysETCCDI','Ann','days'),
                        climdex.r95store=c('r95storeETCCDI','Ann','days'),
                        climdex.r99store=c('r99storeETCCDI','Ann','days'))

  rv <- climdex.names[[climdex.name]]
  return(rv)
}


##-----------------------------------------------------------

create.base.files <- function(var.name,var.units,long.name,
			      gcm,scenario,
                              write.clim.name,
                              time.info,space.info,
                              data.dir,write.dir) {


  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', space.info$lon)
  y.geog <- ncdim_def('lat', 'degrees_north', space.info$lat)
  t.geog <- ncdim_def('time', time.info$atts$units, time.info$dates,
                      unlim=TRUE, calendar=time.info$atts$calendar)
  
  var.geog <- ncvar_def(var.name, units=var.units, dim=list(x.geog, y.geog, t.geog),
                        missval=-32768.)
  print('Create this file')			
  print(paste(write.dir,'/',write.clim.name,sep=''))
  file.nc <- nc_create(paste(write.dir,'/',write.clim.name,sep=''), var.geog)
  
  ##Loop over subsets of the time series
  ##Past file first
  global.atts <- space.info$global.atts
  global.names <- names(global.atts)
  for (g in 1:length(global.names)) {
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])
  }
  print('Finished Global Atts')
  
  ##Time attributes
  ncatt_put(file.nc,varid='time',attname='units',attval=time.info$atts$units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=time.info$atts$calendar)  
  print('Finished time atts')
  
  lon.atts <- space.info$lon.atts
  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts)) {
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])
  }
  print('Finished Lon atts')
  
  lat.atts <- space.info$lat.atts
  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts)) {
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])
  }
  print('Finished lat atts')
  
  ##Climdex Attributes
  ncatt_put(file.nc,varid=var.name,attname='units',attval=var.units)
  ncatt_put(file.nc,varid=var.name,attname='_FillValue',attval=-32768.)
  ncatt_put(file.nc,varid=var.name,attname='standard_name',attval=long.name)
  ncatt_put(file.nc,varid=var.name,attname='long_name',attval=long.name)
  print('Finished Var atts')
  nc_close(file.nc)
}



##**************************************************************************************



##-----------------------------------------------------------------------------
##Create empty degree days files using dimensions from the 800m PRISM adjusted BCCAQ2

make.derived.files <- function(gcm,scenario,
                               tasmax.time.tmp,tasmin.time.tmp,pr.time.tmp,
                               tasmax.space.tmp,tasmin.space.tmp,pr.space.tmp,
                               data.dir,write.dir,tmp.dir) {

  ##-----------------------------------------------------------------------
  ##Degree Days
  
  degree.names <- c('cdd','fdd','gdd','hdd')
  dd.long.names <- c('Cooling Degree Days','Freezing Degree Days',
  		     	  'Growing Degree Days','Heating Degree Days')
  freq <- 'Ann'
  file.split <- strsplit(tasmax.time.tmp,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  out.dir <- paste(tmp.dir,'/degree_days',sep='')
  if (!file.exists(out.dir)) {
    dir.create(out.dir,recursive=TRUE)
  }
  print(tasmax.time.tmp)
  print(tasmax.space.tmp)
  time.info <- time.component(tasmax.time.tmp,freq)
  space.info <- space.component(tasmax.space.tmp)  
  
  for (d in seq_along(degree.names)) {
      degree.name <- degree.names[d]
      long.name <- dd.long.names[d]
      write.clim.name <- paste(degree.name,'_annual_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',time.info$interval,'.nc',sep='')    
      create.base.files(degree.name,'degree_days',long.name,
                        gcm,scenario,    
                        write.clim.name,
                        time.info,space.info,		
                        tmp.dir,out.dir)
  }
  move.back <- paste("rsync -av ",out.dir,' ',write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up.dd <- paste("rm ",tmp.dir,"/degree_days/*nc" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

  clean.up.dd <- paste("rmdir ",tmp.dir,"/degree_days/" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)
  
  ##-----------------------------------------------------------------------
  ##Return Periods
  ##Creating Annual Maximum Values here so any return period year can be calculated
  print('Creating Return Period Annual Maximum base files')
  var.names <- c('pr','tasmax','tasmin')
  var.long.names <- c('Annual Maximum Precipitation',
  		     'Annual Maximum Tasmax',
		     'Annual Minimum Tasmin')
  freq <- 'Ann'                      

  out.dir <- paste(tmp.dir,'/annual_maxima',sep='')
  if (!file.exists(out.dir)) {
    dir.create(out.dir,recursive=TRUE)
  }

  for (v in seq_along(var.names)) {
      var.name <- var.names[v]
      long.name <- var.long.names[v]
    print(var.name)
    time.tmp <- switch(var.name,
                       tasmax=tasmax.time.tmp,
                       tasmin=tasmin.time.tmp,
                       pr=pr.time.tmp)
    print(time.tmp)		        
    time.info <- time.component(time.tmp,freq)
    space.tmp <- switch(var.name,
                        tasmax=tasmax.space.tmp,
                        tasmin=tasmin.space.tmp,
                        pr=pr.space.tmp)
    print(space.tmp)			 
    space.info <- space.component(space.tmp)  
    var.units <- switch(var.name,
                       tasmax='degC',
                       tasmin='degC',
                       pr='mm day-1')
    write.clim.name <- paste(var.name,'_Annual_Max_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',time.info$interval,'.nc',sep='')
    if (var.name=='tasmin') {
        write.clim.name <- paste(var.name,'_Annual_Min_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',time.info$interval,'.nc',sep='')
    }

    create.base.files(var.name,var.units,long.name,
		      gcm,scenario,    
                      write.clim.name,
                      time.info,space.info,		
                      tmp.dir,out.dir)
  }
  move.back <- paste("rsync -av ",out.dir," ",write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up.dd <- paste("rm ",tmp.dir,"/annual_maxima/*nc" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

  clean.up.dd <- paste("rmdir ",tmp.dir,"/annual_maxima/" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

  ##-----------------------------------------------------------------------
  ##Climdex
  print('Creating Climdex Files')
  var.names <- c('climdex.fd','climdex.su','climdex.id',
                 'climdex.tr','climdex.gsl','climdex.txx',
                 'climdex.tnx','climdex.txn','climdex.tnn',
                 'climdex.tn10p','climdex.tx10p','climdex.tn90p',
                 'climdex.tx90p','climdex.wsdi','climdex.csdi',
                 'climdex.dtr','climdex.rx1day','climdex.rx2day',
                 'climdex.rx5day','climdex.sdii','climdex.r10mm',
                 'climdex.r20mm','climdex.cdd','climdex.cwd',
                 'climdex.r95ptot','climdex.r99ptot','climdex.prcptot',
                 'climdex.r95days','climdex.r99days','climdex.r95store',
                 'climdex.r99store')

  out.dir <- paste(tmp.dir,'/climdex',sep='')
  if (!file.exists(out.dir)) {
    dir.create(out.dir,recursive=TRUE)
  }

  for (v in seq_along(var.names)) {
      var.name <- var.names[v]
      climdex.info <- get.climdex.info(var.name)
      climdex.var <- climdex.info[1]
      climdex.calendar <- climdex.info[2]
      climdex.units <- climdex.info[3]

      print(var.name)
      time.tmp <- tasmax.time.tmp
      print(time.tmp)		        
      time.info <- time.component(time.tmp,freq=climdex.calendar)
      space.tmp <- tasmax.space.tmp                 
      print(space.tmp)			 
      space.info <- space.component(space.tmp)  
      var.units <- climdex.units

      write.clim.name <- paste(climdex.var,'_',tolower(climdex.calendar),'_BCCAQ2-PRISM_',gcm,'_',scenario,'_',run,'_',time.info$interval,'.nc',sep='')
      create.base.files(climdex.var,var.units,toupper(climdex.var),
                        gcm,scenario,    
                        write.clim.name,
                        time.info,space.info,		
                        tmp.dir,out.dir)
  }
  move.back <- paste("rsync -av ",out.dir," ",write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up.dd <- paste("rm ",tmp.dir,"/climdex/*nc" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

  clean.up.dd <- paste("rmdir ",tmp.dir,"/climdex/" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)





  ##-----------------------------------------------------------------------
  ##Annual Averages/Sums
  
  var.names <- c('pr','tasmax','tasmin')
  var.long.names <- c('Annual Total Precipitation','Annual Average Maximum Temperature',
  		     	  'Annual Average Minimum Temperature')
  freq <- 'Ann'
  file.split <- strsplit(tasmax.time.tmp,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  out.dir <- paste(tmp.dir,'/annual',sep='')
  if (!file.exists(out.dir)) {
    dir.create(out.dir,recursive=TRUE)
  }
  print(tasmax.time.tmp)
  print(tasmax.space.tmp)
  time.info <- time.component(tasmax.time.tmp,freq)
  space.info <- space.component(tasmax.space.tmp)  
  
  for (v in seq_along(var.names)) {
      var.name <- var.names[v]
      long.name <- var.long.names[v]
      var.units <- switch(var.name,
                          tasmax='degC',
                          tasmin='degC',
                          pr='mm day-1')
      
      write.clim.name <- paste(var.name,'_annual_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',time.info$interval,'.nc',sep='')    
      create.base.files(var.name,var.units,long.name,
                        gcm,scenario,    
                        write.clim.name,
                        time.info,space.info,		
                        tmp.dir,out.dir)
  }
  move.back <- paste("rsync -av ",out.dir,' ',write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up.dd <- paste("rm ",tmp.dir,"/annual/*nc" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

  clean.up.dd <- paste("rmdir ",tmp.dir,"/annual/" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)


##-----------------------------------------------------------------------
  ##Seasonal Averages/Sums
  
  var.names <- c('pr','tasmax','tasmin')
  var.long.names <- c('Seasonal Total Precipitation','Seasonal Average Maximum Temperature',
  		     	  'Seasonal Average Minimum Temperature')
  freq <- 'Seas'
  file.split <- strsplit(tasmax.time.tmp,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  out.dir <- paste(tmp.dir,'/seasonal',sep='')
  if (!file.exists(out.dir)) {
    dir.create(out.dir,recursive=TRUE)
  }
  print(tasmax.time.tmp)
  print(tasmax.space.tmp)
  time.info <- time.component(tasmax.time.tmp,freq)
  space.info <- space.component(tasmax.space.tmp)  
  
  for (v in seq_along(var.names)) {
      var.name <- var.names[v]
      long.name <- var.long.names[v]
      var.units <- switch(var.name,
                          tasmax='degC',
                          tasmin='degC',
                          pr='mm day-1')
      
      write.clim.name <- paste(var.name,'_seasonal_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',time.info$interval,'.nc',sep='')    
      create.base.files(var.name,var.units,long.name,
                        gcm,scenario,    
                        write.clim.name,
                        time.info,space.info,		
                        tmp.dir,out.dir)
  }
  move.back <- paste("rsync -av ",out.dir,' ',write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up.dd <- paste("rm ",tmp.dir,"/seasonal/*nc" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

  clean.up.dd <- paste("rmdir ",tmp.dir,"/seasonal/" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

##-----------------------------------------------------------------------
  ##Monthly Averages/Sums
  
  var.names <- c('pr','tasmax','tasmin')
  var.long.names <- c('Monthly Total Precipitation','Monthly Average Maximum Temperature',
  		     	  'Monthly Average Minimum Temperature')
  freq <- 'Mon'
  file.split <- strsplit(tasmax.time.tmp,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]

  out.dir <- paste(tmp.dir,'/monthly',sep='')
  if (!file.exists(out.dir)) {
    dir.create(out.dir,recursive=TRUE)
  }
  print(tasmax.time.tmp)
  print(tasmax.space.tmp)
  time.info <- time.component(tasmax.time.tmp,freq)
  space.info <- space.component(tasmax.space.tmp)  
  
  for (v in seq_along(var.names)) {
      var.name <- var.names[v]
      long.name <- var.long.names[v]
      var.units <- switch(var.name,
                          tasmax='degC',
                          tasmin='degC',
                          pr='mm day-1')
      
      write.clim.name <- paste(var.name,'_monthly_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',time.info$interval,'.nc',sep='')    
      create.base.files(var.name,var.units,long.name,
                        gcm,scenario,    
                        write.clim.name,
                        time.info,space.info,		
                        tmp.dir,out.dir)
  }
  move.back <- paste("rsync -av ",out.dir,' ',write.dir,sep='')
  print(move.back)
  system(move.back)

  clean.up.dd <- paste("rm ",tmp.dir,"/monthly/*nc" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)

  clean.up.dd <- paste("rmdir ",tmp.dir,"/monthly/" ,sep='')
  print(clean.up.dd)      
  system(clean.up.dd)




}


##****************************************************************************************
##****************************************************************************************

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}
tmp.dir <- tmpdir

##tmp.dir <- '/tmpdir'
gcm <- 'CanESM2'
scenario <- 'rcp85'
data.dir <- '/scratch/ssobie/prism/'
write.dir <- paste0('/scratch/ssobie/prism/',gcm,'/',scenario)

##Move data to local storage for better I/O
if (!file.exists(tmp.dir))
  dir.create(tmp.dir,recursive=TRUE)    

##Transfer template files for derived file creation
print('Transfer Files')
tasmax.time.file <- list.files(paste0(data.dir,gcm,'/lat_split/'),pattern=paste0('tasmax_gcm_prism'))[1]
tasmax.time.tmp <- paste0(tmp.dir,'/',gcm,'/',tasmax.time.file)
move.to <- paste("rsync -av ",paste0(data.dir,gcm,'/lat_split/',tasmax.time.file),' ',tmp.dir,'/',gcm,'/',sep='')
print(move.to)
system(move.to)

tasmin.time.file <- list.files(paste0(data.dir,gcm,'/lat_split/'),pattern=paste0('tasmin_gcm_prism'))[1]
tasmin.time.tmp <- paste0(tmp.dir,'/',gcm,'/',tasmin.time.file)
move.to <- paste("rsync -av ",paste0(data.dir,gcm,'/lat_split/',tasmin.time.file),' ',tmp.dir,'/',gcm,'/',sep='')
print(move.to)
system(move.to)

pr.time.file <- list.files(paste0(data.dir,gcm,'/lat_split/'),pattern=paste0('pr_gcm_prism'))[1]
pr.time.tmp <- paste0(tmp.dir,'/',gcm,'/',pr.time.file)
move.to <- paste("rsync -av ",paste0(data.dir,gcm,'/lat_split/',pr.time.file),' ',tmp.dir,'/',gcm,'/',sep='')
print(move.to)
system(move.to)

tasmax.space.file <- list.files(paste0(data.dir,gcm,'/daily_prism/'),pattern=paste0('tasmax_gcm_prism'))[1]
tasmax.space.tmp <- paste0(tmp.dir,'/',gcm,'/',tasmax.space.file)
move.to <- paste("rsync -av ",paste0(data.dir,gcm,'/daily_prism/',tasmax.space.file),' ',tmp.dir,'/',gcm,'/',sep='')
print(move.to)
system(move.to)

tasmin.space.file <- list.files(paste0(data.dir,gcm,'/daily_prism/'),pattern=paste0('tasmin_gcm_prism'))[1]
tasmin.space.tmp <- paste0(tmp.dir,'/',gcm,'/',tasmin.space.file)
move.to <- paste("rsync -av ",paste0(data.dir,gcm,'/daily_prism/',tasmin.space.file),' ',tmp.dir,'/',gcm,'/',sep='')
print(move.to)
system(move.to)

pr.space.file <- list.files(paste0(data.dir,gcm,'/daily_prism/'),pattern=paste0('pr_gcm_prism'))[1]
pr.space.tmp <- paste0(tmp.dir,'/',gcm,'/',pr.space.file)
move.to <- paste("rsync -av ",paste0(data.dir,gcm,'/daily_prism/',pr.space.file),' ',tmp.dir,'/',gcm,'/',sep='')
print(move.to)
system(move.to)

make.derived.files(gcm,scenario,
                   tasmax.time.tmp,tasmin.time.tmp,pr.time.tmp,
                   tasmax.space.tmp,tasmin.space.tmp,pr.space.tmp,
                   data.dir,write.dir,tmp.dir)

clean.up <- paste("rm ",tmp.dir,'/',gcm,"/*nc" ,sep='')
print(clean.up)      
system(clean.up)
