##******************************************************************************
##******************************************************************************

library(ncdf4)
library(PCICt)
library(fields)


##******************************************************************************
################################################################################

daily.prism.scale <- function(var.name,prism.var,gcm,year,prism.file,base.dir,tmp.dir) {

  interval <- year

  print(paste('Daily PRISM: ',gcm,', ',var.name,', ',interval,sep=''))
  gcm.dir <- paste(base.dir,gcm,'/',sep='')
  bccaq.files <- list.files(path=paste0(gcm.dir,'interpolated/'),
			    pattern=paste(var.name,'_anoms_interp_BCCAQ2_',sep=''))
  bccaq.file <- bccaq.files[grep(interval,bccaq.files)]
  bccaq.tmp <- paste0(tmp.dir,bccaq.file)

  move.to <- paste0('rsync -av ',gcm.dir,'interpolated/',bccaq.file,' ',tmp.dir)
  print(move.to)
  system(move.to)

  adjusted.file <- gsub(pattern='_anoms_interp_BCCAQ2_',replacement='_gcm_prism_BCCAQ2_',bccaq.tmp)
  print('Adjusted File')		
  print(adjusted.file)  		
  file.copy(from=bccaq.tmp,to=adjusted.file,overwrite=T)  

  ##PRISM climatologies
  print('PRISM files')
  print(prism.file)
  pnc <- nc_open(prism.file)
  prism.clim <- ncvar_get(pnc,prism.var)
  nc_close(pnc)

  print(paste0('Open file: ',bccaq.tmp))
  bnc <- nc_open(bccaq.tmp)
  print(paste0('Open file: ',adjusted.file))  
  anc <- nc_open(adjusted.file,write=TRUE)  
  print('Opened both files')
  
  time.atts <- ncatt_get(bnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(bnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400
  monthly.fac <- as.factor(format(var.dates,'%m'))

  data.atts <- ncatt_get(anc,var.name)
  print(data.atts)
  nlon <- bnc$dim$lon$len
  nlat <- bnc$dim$lat$len
  addon <- matrix(NA,nrow=nlon,ncol=9)
  
  var.data <- ncvar_get(bnc,var.name)
  nc_close(bnc)
  var.adjust <- var.data*0  
  fac.sub <- monthly.fac
  for(mn in 1:12) {
    print(mn)
    prism.mean <- prism.clim[,,mn] ##cbind(addon,prism.clim[,,mn])
    var.ix <- which(fac.sub==sprintf('%02d',mn))
    mlen <- length(var.ix)
    for (i in 1:mlen) {
      ix <- var.ix[i]
      var.sub <- var.adjust[,,ix]
      if (var.name=='pr') {
        var.sub <- var.data[,,ix]*prism.mean
      }
      if (grepl('tas',var.name)) {
        var.sub <- var.data[,,ix] + prism.mean
      }
      var.adjust[,,ix] <- var.sub
    }##Loop over indices    
  }##Loop over Months
  rm(var.sub)
  rm(var.data)
  print(range(var.adjust,na.rm=T))
  ##var.adjust[var.adjust>800] <- 777
  ##data.adjust <- (var.adjust - data.atts$add_offset)/data.atts$scale_factor
  ncvar_put(anc,varid=var.name,vals=var.adjust)
  rm(var.adjust)
  nc_close(anc)

  move.back <- paste0('rsync -av ',adjusted.file,' ',gcm.dir,'daily_prism/')
  print(move.back)
  system(move.back)

  clean.up <- paste0('rm ',bccaq.tmp)
  print(clean.up)
  system(clean.up)

  clean.up <- paste0('rm ',adjusted.file)
  print(clean.up)
  system(clean.up)

  gc()
}
################################################################################
##******************************************************************************

run.adjust <- function() {

  ##Requires 1981-2010 (or equivalent base period) in the /baseline directory to create anomalies from the full period
  ## (usually 1950-2000 and 2001-2100). Both of these are created using extract.bccaq.gcm.r
  ##Also need the PRISM climatologies (also using extract.bccaq.gcm.r).

  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
  var.name <- varname

  tmp.dir <- tmpdir
  if (!file.exists(tmp.dir)) {
     dir.create(tmp.dir,recursive=T)
  }

  base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
  years <- seq(1951,1955,by=1)
    
  prism.var <- switch(var.name,
                      pr='pr',
                      tasmax='tmax',
                      tasmin='tmin')
  
  prism.file <- paste(base.dir,'PRISM/',prism.var,'_day_PRISM_observation_1981-2010.nc',sep='')
  for (year in years) {
    daily.prism.scale(var.name,prism.var,gcm,year,prism.file,base.dir,tmp.dir)
  }   

  clean.up <- paste("rm ",tmp.dir,"/*nc" ,sep='')
  print(clean.up)
  system(clean.up)

}

run.adjust()



