##******************************************************************************
##******************************************************************************
library(ncdf4)
library(PCICt)
library(fields)

##******************************************************************************
################################################################################
##Calculate BCCAQ anomalies

create.anoms <- function(var.name,file,var.mon) {

  ##For mean subset 1981-2010
  anoms.file <- gsub(pattern='_day_',replacement='_anoms_',file)
  print('Anoms file')
  print(anoms.file)
  file.copy(from=file,to=anoms.file,overwrite=T)
  Sys.sleep(5)
  if (!file.exists(anoms.file))
    browser()
  
  nc <- nc_open(file)
  anc <- nc_open(anoms.file,write=TRUE)
  lat <- ncvar_get(nc,'lat') ##n.lat = 145
  l.seq <- seq(1,145,by=5)
  lon <- ncvar_get(nc,'lon')
  n.lon <- length(lon)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  n.time <- length(time.values)
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400
  monthly.fac <- as.factor(format(var.dates,'%m'))

  data.atts <- ncatt_get(nc,var.name)
  for (ltx in l.seq) { ##145/5=29
    print(paste0('Subset: ',ltx,' in 145'))
    var.data <- ncvar_get(nc,var.name,start=c(1,ltx,1),count=c(-1,5,-1))
    
    lt.ix <- ltx:(ltx+4)
    ##Load full time series and take anomalies from this
    var.anoms <- var.data*0
    for(mn in 1:12) {
      ##print(mn)
      var.mean <- var.mon[mn,,lt.ix]
      var.ix <- which(monthly.fac==sprintf('%02d',mn))
      mlen <- length(var.ix)
      for (i in 1:mlen) {
        ix <- var.ix[i]
        if (var.name=='pr')
          var.anoms[,,ix] <- var.data[,,ix]/var.mean
        if (grepl('tas',var.name))
          var.anoms[,,ix] <- var.data[,,ix] - var.mean

      }
    }

    ##Fix the edges for interpolations
    if (1==1) {
    ncol <- dim(var.data)[2]
    for (l in 1:3) { ##Repeat 3 times to add buffer
      for (j in 1:(ncol-1)) {
        ix <- is.na(var.anoms[,j,1])
        var.anoms[ix,j,] <- var.anoms[ix,(j+1),]    
      }
      for (j in ncol:2) {
        ix <- is.na(var.anoms[,j,1])
        var.anoms[ix,j,] <- var.anoms[ix,(j-1),]    
      }

      nrow <- dim(var.data)[1]
      for (k in 1:(nrow-1)) {
        kx <- is.na(var.anoms[k,,1])
        var.anoms[k,kx,] <- var.anoms[(k+1),kx,]    
      }
      for (k in nrow:2) {
        kx <- is.na(var.anoms[k,,1])
        var.anoms[k,kx,] <- var.anoms[(k-1),kx,]    
      }
    }
    }
    ##if (var.name=='pr')
    ##  var.anoms[is.na(var.anoms)] <- 0
    ##data.adjust <- (var.anoms - data.atts$add_offset)/data.atts$scale_factor

    ncvar_put(anc,varid=var.name,vals=var.anoms,start=c(1,ltx,1),count=c(n.lon,5,n.time))
  
  }
  nc_close(nc)
  nc_close(anc)
  gc()
  return(anoms.file)
}

bccaq.anomalies <- function(var.name,gcm,scenario,base.dir,tmp.dir) {
  print(paste('BCCAQ Anomalies: ',gcm,', ',var.name,sep=''))

  base.files <- list.files(path=paste(base.dir,'baseline/',gcm,sep=''),pattern=paste(var.name,'_day_',sep=''))
  base.file <- base.files[grep('rcp85',base.files)] 

  move.to <- paste("rsync -av ",base.dir,"baseline/",gcm,"/",base.file," ",tmp.dir,sep='')
  print(move.to)
  system(move.to)

  base.file <- paste0(tmp.dir,'/',base.file)
  print(base.file)

  gnc <- nc_open(base.file)
  time.atts <- ncatt_get(gnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(gnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400
  monthly.fac <- as.factor(format(var.dates,'%m'))
  monthly.ts.fac <- as.factor(format(var.dates,'%Y-%m'))
  mon.facs <- as.factor(format(as.Date(paste(levels(monthly.ts.fac),'-01',sep='')),'%m'))
  var.data <- ncvar_get(gnc,var.name)
  var.mon <- apply(var.data,c(1,2),function(x,fac){tapply(x,fac,mean,na.rm=T)},monthly.fac)
  if (var.name=='pr') {
    var.data[var.data <=0] <- NA    
    var.test <- apply(var.data,c(1,2),function(x,fac){tapply(x,fac,sum,na.rm=T)},monthly.ts.fac)
    var.mon <-  apply(var.test,c(2,3),function(x,fac){tapply(x,fac,mean,na.rm=T)},mon.facs)
  }
  nc_close(gnc)
  ##------------------------------------  
  gcm.dir <- paste(base.dir,gcm,sep='')
  var.files <- list.files(path=gcm.dir,pattern=paste(var.name,'_day_BCCAQ2',sep=''))
  scen.files <- var.files[grep(scenario,var.files)]
  past.file <- scen.files[grep('1951-2000',scen.files)]
  proj.file <- scen.files[grep('2001-2100',scen.files)]

  move.to <- paste("rsync -av ",gcm.dir,"/",past.file," ",tmp.dir,sep='')
  print(move.to)
  system(move.to)
  
  move.to <- paste("rsync -av ",gcm.dir,"/",proj.file," ",tmp.dir,sep='')
  print(move.to)
  system(move.to)
  
  past.file <- paste0(tmp.dir,'/',past.file)
  proj.file <- paste0(tmp.dir,'/',proj.file)

  anoms.past <- create.anoms(var.name,past.file,var.mon)
  move.back <- paste("rsync -av ",anoms.past," ",gcm.dir,sep='')
  print(move.back)
  system(move.back)

  anoms.proj <- create.anoms(var.name,proj.file,var.mon)
  move.back <- paste("rsync -av ",anoms.proj," ",gcm.dir,sep='')
  print(move.back)
  system(move.back)

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

##gcm <- 'inmcm4'
##scenario <- 'rcp85'
##varname <- 'tasmin'

  tmp.dir <- paste0('/local_temp/ssobie/anomalies/',gcm,'/',varname,'/')
  if (!file.exists(tmp.dir)) {
     dir.create(tmp.dir,recursive=T)
  }


  base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
  grid.file <- '/storage/home/ssobie/grid_files/bc.prism.grid.txt'

  bccaq.anomalies(varname,gcm,scenario,base.dir,tmp.dir)

  clean.up <- paste("rm ",tmp.dir,"/*nc" ,sep='')
  print(clean.up)
  system(clean.up)

}

ptm <- proc.time()
run.adjust()
print('Elapsed time:')
print(proc.time()-ptm)



