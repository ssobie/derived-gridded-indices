##Script to convert PCSD ASCII grid to a netcdf format

library(ncdf4)
library(PCICt)

##**********************************************************************
##----------------------------------------------------------------------

generate_monthly_time_series <- function(startdate,enddate) {
  months <- seq(from=as.Date(startdate),by='month',to=as.Date(enddate))
}

##----------------------------------------------------------------------
##Global Attributes
get_global_atts <- function() {
  global.atts <- list(institution="Pacific Climate Impacts Consortium (PCIC), Victoria, BC, www.pacificclimate.org",
                   contact="Pacific Climate Impacts Consortium",
                   Conventions="CF-1.4",
                   institute_id ="PCIC",
                   domain='Canada',
                   creation_date=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%Z'),
                   frequency="month",
                   product="Gridded Provincial Climate Dataset",
                   modeling_realm="atmos",
                   project_id='PCDS',
                   references=paste0("Pacific Climate Impacts Consortium, ",format(Sys.time(),'%Y'),": Provincial climate dataset. Pacific Climate Impacts Consortium [Available online at https://www.pacificclimate.org/data/bc-station-data.]"),
                   title = "Gridded Provincial Climate dataset")
  return(global.atts)
}

##----------------------------------------------------------------------
get_standard_atts <- function(var.name) {
  lon.atts <- list(standard_name="longitude",
                   long_name = "longitude",
                   units = "degrees_east",
                   axis = "X")
  lat.atts <- list(standard_name="latitude",
                   long_name = "latitude",
                   units = "degrees_north",
                   axis = "Y")
  pr.atts <- list(standard_name = "precipitation",
                  long_name = "Precipitation",
                  missing_value = -9999.0,
                  cell_methods = "time: sum",
                  units = "%")
  tasmax.atts <- list(standard_name = "air_temperature",
                      long_name = "Daily Maximum Near-Surface Air Temperature",
                      units = "degC",
                      missing_value = -9999.0,
                      cell_methods = "time: maximum")
  tasmin.atts <- list(standard_name = "air_temperature",
                      long_name = "Daily Minimum Near-Surface Air Temperature",
                      units = "degC",
                      missing_value = -9999.0,
                      cell_methods = "time: minimum")
  var.atts <- switch(var.name,
                     pr=pr.atts,
                     tasmax=tasmax.atts,
                     tasmin=tasmin.atts)

  rv <- list(lon=lon.atts,
             lat=lat.atts,
             var=var.atts)
  return(rv)
}

##----------------------------------------------------------------------


make_pcds_netcdf <- function(var.name,pcds.array,                            
                             lon,lat,
                             startdate,enddate,
                             file.dir,file.name) {
  ##Time 
  time.calendar <- 'gregorian'
  time.units <- paste0('days since ',startdate)
  time.series <- generate_monthly_time_series(startdate,enddate)
  time.values <- as.numeric(time.series - as.Date(startdate))
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  t.geog <- ncdim_def('time', time.units, time.values,
                      unlim=TRUE, calendar=time.calendar)

  var.units <- list(pr='%',
                    tasmax='degC',
                    tasmin='degC')

  ##Netcdf dimensions
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  var.geog <- ncvar_def(varname, units=var.units[[varname]], dim=list(x.geog, y.geog, t.geog),
                        missval=-9999.0)
  pcds.nc <- nc_create(paste(file.dir,file.name,sep=''), var.geog)
  
  ##Time Attributes
  ncatt_put(pcds.nc,'time','standard_name','Time')
  ncatt_put(pcds.nc,'time','long_name','Time')

  ##Standard Attributes
  atts <- get_standard_atts(varname)
  print('Lon names')
  lon.names <- names(atts$lon)
  for (j in 1:length(atts$lon))
    ncatt_put(pcds.nc,varid='lon',attname=lon.names[j],attval=atts$lon[[j]])
  print('Lat names')
  lat.names <- names(atts$lat)
  for (j in 1:length(atts$lat))
    ncatt_put(pcds.nc,varid='lat',attname=lat.names[j],attval=atts$lat[[j]])
  print('Var names')
  var.names <- names(atts$var)
  for (j in 1:length(atts$var))
    ncatt_put(pcds.nc,varid=varname,attname=var.names[j],attval=atts$var[[j]])

  ##Global Attributes
  global.atts <- get_global_atts()
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(pcds.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ncvar_put(pcds.nc,'lon',lon)
  ncvar_put(pcds.nc,'lat',lat)
  ncvar_put(pcds.nc,varname,pcds.array)
    
  nc_close(pcds.nc)
  return(paste0(file.dir,file.name))
}


##----------------------------------------------------------------------

read_gridded_pcds <- function(file.dir,file.name) {
  file <- paste0(file.dir,file.name)
  ##input <- as.matrix(read.table(gzfile(file),skip=6))
  ##header <- read.table(gzfile(file),nrows=6)
  input <- as.matrix(read.table(file,skip=6))
  header <- read.table(file,nrows=6)
##browser()
  rv <- list(data=input,header=header)
  return(rv)
}

##----------------------------------------------------------------------   

rotate <- function(x) t(apply(x, 2, rev))

assemble_pcds_layers <- function(varname,lon,lat,startdate,enddate,read.dir) {

  var.short <- switch(varname,pr='ppt',tasmax='tx',tasmin='tn')

  monthly.series <- generate_monthly_time_series(startdate,enddate) 
  pcds.array <- array(0,c(length(lon),length(lat),length(monthly.series)))  

  for (i in seq_along(monthly.series)) {
     month <- as.numeric(format(monthly.series[i],'%m'))
     year <- format(monthly.series[i],'%Y')
     file.dir <- paste0(read.dir,month,'_',var.short,'_anoms/')
     file.name <- paste0('ahccd_D_ann_annmean_',year,'_NA_se.dat')
     ##file.name <- paste0('ahccd_D_ann_annmean_',year,'_NA.dat.gz')
     pcds.read <- read_gridded_pcds(file.dir,file.name)
     pcds.array[,,i] <- rotate(pcds.read$data)
  }
  if (varname=='pr') 
     pcds.array <- pcds.array*100

  return(pcds.array)
}

##----------------------------------------------------------------------
##**********************************************************************


var.list <- c('pr','tasmax','tasmin')
##Spatial Coordinates
lon <- seq(-145.0,-105.0,0.5)
lat <- seq(  48.0,  61.0,0.5)

startdate <- '1900-01-01'
enddate <- '2017-12-31'

for (varname in var.list) {

file.dir <- '/storage/data/projects/rci/data/assessments/bc/pcds/'
file.name <- paste0(varname,'_PCDS_',format(as.Date(startdate),'%Y'),'-',       
                                     format(as.Date(enddate),'%Y'),'.SE.nc')
read.dir <- '/storage/data/projects/crmp/bc_trends_indicators/'

pcds.series <- assemble_pcds_layers(varname,
                                    lon,lat,
                                    startdate,enddate,
                                    read.dir)

pcds.file <- make_pcds_netcdf(varname,pcds.series,                            
                              lon,lat,
                              startdate,enddate,
                              file.dir,file.name)

}
