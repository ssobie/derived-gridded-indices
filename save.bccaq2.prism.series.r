##Contains functions to compute Building Code values from the BC
##Building Code and the American ASHRAE Handbook

library(ncdf4)
library(PCICt)

get.time.series <- function(nc) {

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  values <- ncvar_get(nc,'time')
  series <- format(origin + (values)*86400,'%Y-%m-%d')
  return(series)
}

read.data.subset <- function(data,dates,interval) {

  bnds <- strsplit(interval,'-')[[1]]
  st <- head(grep(bnds[1],dates),1)
  en <- tail(grep(bnds[2],dates),1)
  if (length(en)==0) {
    en <- length(dates)
  }
  data.subset <- data[st:en]
  dates.subset <- dates[st:en]
  rv <- data.subset
  return(rv)
}

get.bccaq.data <- function(var.name,gcm,scenario,lats,data.dir,lon.c,lat.c) {

   var.files <- list.files(path=data.dir,pattern=var.name,full.name=TRUE)
   data.file <- var.files[grep(lats,var.files)]
   print(data.file)
   nc <- nc_open(data.file)
 
   lon <- ncvar_get(nc,'lon')
   lon <- ((lon + 180) %% 360) - 180
   lat <- ncvar_get(nc,'lat')
  
   lon.ix <- which.min(abs(lon.c-lon))
   lat.ix <- which.min(abs(lat.c-lat))

   data.subset <- ncvar_get(nc,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
   time.series <- get.time.series(nc)
   time.subset <- as.Date(time.series)

   nc_close(nc)

   rv <- list(data=data.subset,time=time.subset)
   return(rv)
}


##BCCAQ Data
gather.bccaq.data <- function(gcm,lon.c,lat.c,lats,scenario) {

  bccaq.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',gcm,'/lat_split/')

  tasmax.data <- get.bccaq.data('tasmax',gcm,scenario,lats,bccaq.dir,lon.c,lat.c)
  tasmin.data <- get.bccaq.data('tasmin',gcm,scenario,lats,bccaq.dir,lon.c,lat.c)
  tas.data <- tasmax.data
  tas.data$data <- (tasmax.data$data+tasmin.data$data)/2
  pr.data <- get.bccaq.data('pr',gcm,scenario,lats,bccaq.dir,lon.c,lat.c)

  rv <- list(tasmax=tasmax.data,        
             tasmin=tasmin.data,
             tas=tas.data,
             pr=pr.data)
  return(rv)
}


##************************************************************************
##Load GCM BasedData


gcm <- 'CSIRO-Mk3-6-0'

scenario <- 'rcp85'



##Lions Gate Hospital
lon.c <- -123.0685
lat.c <- 49.3210
lats <- '49.3-49.4'

##BCCAQ Data
bccaq.data <- gather.bccaq.data(gcm,lon.c,lat.c,lats,scenario)

save.dir <- '/storage/data/climate/downscale/CMIP5/building_code/data_files/'
save(bccaq.data,file=paste0(save.dir,'lions_gate_hospital_bccaq2_',gcm,'_rcp85_pr_tx_tn.RData'))

##Bella Bella
lon.c <- -128.1434
lat.c <- 52.1615
lats <- '52.1-52.2'

##BCCAQ Data
bccaq.data <- gather.bccaq.data(gcm,lon.c,lat.c,lats,scenario)

save.dir <- '/storage/data/climate/downscale/CMIP5/building_code/data_files/'
save(bccaq.data,file=paste0(save.dir,'bella_bella_hospital_bccaq2_',gcm,'_rcp85_pr_tx_tn.RData'))