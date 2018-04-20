##******************************************************************************
##******************************************************************************

ftm <- proc.time()

library(ncdf4)
library(PCICt)
library(fields)


##******************************************************************************
################################################################################

daily.prism.scale <- function(var.name,prism.var,gcm,year,prism.clim,base.dir,read.dir,tmp.dir) {

  interval <- year

  print(paste('Daily PRISM: ',gcm,', ',var.name,', ',interval,sep=''))
  gcm.dir <- paste(base.dir,gcm,'/',sep='')
  bccaq.files <- list.files(path=paste0(read.dir,'interpolated/'),pattern=paste(var.name,'_anoms_interp_BCCAQ2_',gcm,sep=''))
  bccaq.file <- bccaq.files[grep(interval,bccaq.files)]
  bccaq.tmp <- paste0(tmp.dir,bccaq.file)  
  file.copy(from=paste0(read.dir,'interpolated/',bccaq.file),to=tmp.dir)

  adjusted.file <- gsub(pattern='_anoms_interp_BCCAQ2_',replacement='_gcm_prism_BCCAQ2_',bccaq.tmp)
  print('Adjusted File')		
  print(adjusted.file)  		
  file.copy(from=bccaq.tmp,to=adjusted.file,overwrite=T)    

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
  ntime <- bnc$dim$time$len

  ##Split into loop to handle memory issues effectively
  l.st <- seq(1,1440,by=40)
  l.en <- c(seq(40,1400,by=40),1441)
  l.cnt <- l.en-l.st+1
  stm <- proc.time()    
  for (j in seq_along(l.st)) {

    var.data <- ncvar_get(bnc,var.name,start=c(1,l.st[j],1),count=c(-1,l.cnt[j],-1))
    var.adjust <- var.data*0  
    fac.sub <- monthly.fac
    for(mn in 1:12) {
      print(mn)
      prism.mean <- prism.clim[,(l.st[j]:l.en[j]),mn] 
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
        rm(var.sub)
      }##Loop over indices    
    }##Loop over Months
    rm(var.data)
    ncvar_put(anc,varid=var.name,vals=var.adjust,start=c(1,l.st[j],1),count=c(nlon,l.cnt[j],ntime))
    rm(var.adjust)
  }
  print('Split loop time')
  print(proc.time() - stm)
  nc_close(bnc)
  nc_close(anc)
  
  print('Move back')
  print(adjusted.file)
  print(paste0(read.dir,'daily_prism/'))
  file.copy(from=adjusted.file,to=paste0(read.dir,'daily_prism/'),overwrite=T)            

  file.remove(bccaq.tmp)
  file.remove(adjusted.file)

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

###gcm <- 'HadGEM2-ES'
###scenario <- 'rcp85'
###varname <- 'tasmax'
###tmpdir <- '/local_temp/ssobie/test/'
var.name <- varname
print(varname)

  tmp.dir <- paste0(tmpdir,gcm,'/',varname,'/')
  if (!file.exists(tmp.dir)) {
     dir.create(tmp.dir,recursive=T)
  }

  base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
  read.dir <- paste0('/storage/data/climate/downscale/CMIP5_delivery/',gcm,'/') 
  ##read.dir <- paste0(base.dir,gcm,'/')

  years <- seq(1951,2100,by=1)
  prism.var <- switch(var.name,
                      pr='pr',
                      tasmax='tmax',
                      tasmin='tmin')
  
  prism.file <- paste(base.dir,'PRISM/',prism.var,'_day_PRISM_observation_1981-2010.nc',sep='')
  ##PRISM climatologies
  print('PRISM files')
  print(prism.file)
  pnc <- nc_open(prism.file)
  prism.clim <- ncvar_get(pnc,prism.var)
  nc_close(pnc)

  for (year in years) {
    daily.prism.scale(var.name,prism.var,gcm,year,prism.clim,base.dir,read.dir,tmp.dir)
  }   
  ##file.remove(paste0(tmp.dir,'/*.nc'))
  ##clean.up <- paste("rm ",tmp.dir,"/*nc" ,sep='')
  ##print(clean.up)
  ##system(clean.up)

}



run.adjust()

print('Elapsed time')
print(proc.time()-ftm)

