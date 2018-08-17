##Script to plot regional time series of precipitation and temperature data
##Should add the boundaries of the 8 resource regions to the plot as well.

library(sp)
library(raster)
library(rgdal)

library(ncdf4)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)

library(PCICt)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##-------------------------------------------------------------------------
##-------------------------------------------------------------------------
##Shape and boundary functions

get.region.shape <- function(region,shape.dir) {
  region.shp <- readOGR(shape.dir, region, stringsAsFactors=F, verbose=F)
  return(region.shp)
}

clip.shape.using.raster <- function(file.name,var.name,gcm,region.shp) {

   file.brick <- brick(file.name)
   cells <- cellFromPolygon(file.brick,region.shp,weights=FALSE)
   rsub <- rasterFromCells(file.brick,cells[[1]],values=FALSE)
   
   mon.series <- foreach(i=1:dim(file.brick)[3],.packages='raster',.combine=rbind,.inorder=TRUE) %dopar% {
                         cellStats(crop(file.brick[[i]],rsub),mean,na.rm=T) ##getValues(crop(file.brick[[i]],rsub))
   }   

   rv <- mon.series
   return(rv)
}


##-------------------------------------------------------------------------
##-------------------------------------------------------------------------
monthly.series <- function(region,var.name,gcm.list,scenario) {

  
  ##Directories
  tmp.dir <- '/local_temp/ssobie/monthly/'
  shape.dir <- paste('/storage/data/projects/rci/data/assessments/shapefiles/',region,'/',sep='')   
  interval <- '1951-2100'
  
  ##-------------------------------------------------------------------------
  
  region.shp <- spTransform(get.region.shape(region,shape.dir),CRS("+init=epsg:4326")) ##Keep this projection to extract the data from lat/lon

  ##Prep for the ensemble files
  yrs <- 150
  month.all <- matrix(NA,nrow=yrs,ncol=length(gcm.list))
  sub.names <- toupper(month.abb)
  rv.all <- vector(mode='list',length=12)
  for (m in 1:12)
    rv.all[[m]] <- month.all
  
  for (g in seq_along(gcm.list)) {
    gcm <- gcm.list[g]
    read.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',gcm,'/rcp85/monthly/')
    print(gcm)
    mon.file <- list.files(path=read.dir,pattern=var.name)
    print('Copying')
    file.copy(from=paste0(read.dir,mon.file),to=tmp.dir)
    print('Done copying')
    tmp.file <- paste0(tmp.dir,mon.file)     ##paste0(read.dir,mon.file) ## 

    ##Data clipped by region shapefile
    rtm <- proc.time()
    mon.data <- clip.shape.using.raster(tmp.file,var.name,gcm,region.shp)
    mon.series <- mon.data
    print('Regional clip time')
    print(proc.time()-rtm)

    nc <- nc_open(tmp.file)
    ts <- netcdf.calendar(nc)
    mon.ts <- format(ts,'%m')
    nc_close(nc)

    for (s in seq_along(sub.names)) {
      mon.ix <- mon.ts %in% sprintf('%02d',s)
      mon.month <- mon.series[mon.ix]
      if (length(mon.month) > 150) 
         browser()
      if (length(mon.month) < 150) {
         mon.month <- c(mon.month,NA) 
      }
      rv.all[[s]][,g] <- mon.month
    }

    file.remove(tmp.file)
  }##GCM Loop

  return(rv.all)
}

###***********************************************************************************

  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }

##region <- 'van_coastal_health'
##area <- 'van_coastal_health'
##varname <- 'tasmax'
##tmpdir <- '/local_temp/ssobie/monthly/'

gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

tmp.dir <- tmpdir
if (!file.exists(tmp.dir)) 
  dir.create(tmp.dir)
var.name <- varname

save.dir <- paste0('/storage/data/projects/rci/data/assessments/',area,'/monthly_data_files/')
rcp85.file <- paste(save.dir,var.name,'_',region,'_rcp85_monthly_time_series.RData',sep='')
data.rcp85 <- monthly.series(region,var.name,
                              gcm.list,'rcp85')
save(data.rcp85,file=rcp85.file)

