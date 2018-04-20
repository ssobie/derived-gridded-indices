##Script to calculate and write the standard set of derived variables
##for the 800m data

ptm <- proc.time()
cedar <- FALSE

if (cedar) {
   source('/home/ssobie/assessments/new.netcdf.calendar.R')
   source('/home/ssobie/assessments/cedar.derived.variable.support.r')
} else {
   source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
   source('/home/ssobie/code/repos/derived_gridded_indices/cedar.derived.variable.support.r')
}

library(ncdf4)
library(PCICt)
library(climdex.pcic)

library(foreach)

get.seasonal.fac <- function(seas.dates) {
   seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
   years <- as.numeric(format(seas.dates,'%Y'))
   uni.yrs <- unique(years)
   months <- as.numeric(format(seas.dates,'%m'))
   dec.ix <- grep(12,months)
   years[dec.ix] <- years[dec.ix] + 1
   dec.fix <- years %in% uni.yrs
   yearly.fac <- as.factor(years[dec.fix])
   monthly.fac <- as.factor(format(seas.dates[dec.fix],'%m'))      
   seasonal.fac <- factor(seasons[monthly.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
   avg.fac <- list(yearly.fac,seasonal.fac)       

  return(list(fac=avg.fac,fix=dec.fix))
}

##----------------------------------------------------------------------------------------------
annual.averages.for.model <- function(ann.name,ann.ncs,lat.ix,n.lon,yearly.fac,flag,
                                      ann.list) {                                            
   ##Variables
       flen <- sum(!flag)
       ann.fx <- ann.fxns[[ann.name]]
       sub.list <- ann.list[!flag]
       rm(ann.list)
       
       ann.avg.values <- foreach(
                         data=sub.list,
                         .export=c('yearly.fac','ann.fx')
                         ) %do% {
                              ann.avg.values <- ann.fx(data,yearly.fac)
                         }
      rm(sub.list)                         
      ncol <- length(ann.avg.values[[1]])
      ann.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(ann.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
      rm(ann.avg.values)
      ann.avg.matrix[!flag] <- sub.matrix
      rm(sub.matrix)
      ncvar_put(ann.ncs,varid=ann.name,vals=ann.avg.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
      rm(ann.avg.matrix)                
      gc()
}

##----------------------------------------------------------------------------------------------
seasonal.averages.for.model <- function(seas.name,seas.ncs,lat.ix,n.lon,seasonal.fac,flag,
                                        seas.list) {

       ##Roll over the December months to compute proper seasons
       flen <- sum(!flag)
       seas.nc <- seas.ncs
       seas.fx <- seas.fxns[[seas.name]]
       avg.fac <- seasonal.fac$fac
       dec.fix <- seasonal.fac$fix
       seas.corr.list <- lapply(seas.list,function(x,y){x[y]},dec.fix)
       rm(seas.list)
       sub.list <- seas.corr.list[!flag]
       rm(seas.corr.list)
       seas.avg.values <- foreach(
                         data=sub.list,
                         .export=c('avg.fac','seas.fx')
                         ) %do% {
                              seas.avg.values <- as.vector(t(seas.fx(data,avg.fac)))
                         }
      rm(sub.list)                         		 
      ncol <- length(seas.avg.values[[1]])
      seas.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
      sub.matrix <- matrix(unlist(seas.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
      rm(seas.avg.values)
      seas.avg.matrix[!flag,] <- sub.matrix
      rm(sub.matrix)
      ncvar_put(seas.ncs,varid=seas.name,vals=seas.avg.matrix,
                start=c(1,lat.ix,1),count=c(-1,1,-1))
      rm(seas.avg.matrix)          
      gc()
}

##----------------------------------------------------------------------------------------------
monthly.averages.for.model <- function(mon.name,mon.ncs,lat.ix,n.lon,monthly.fac,flag,
                                       mon.list) {
      
   ##Variables
   flen <- sum(!flag)
   mon.fx <- mon.fxns[[mon.name]]
   sub.list <- mon.list[!flag]
   rm(mon.list)

   mon.avg.values <- foreach(
                             data=sub.list,
                             .export=c('monthly.fac','mon.fx')
                             ) %do% {
                                 mon.avg.values <- mon.fx(data,monthly.fac)
                             }
   rm(sub.list)                         
			 
   ncol <- length(mon.avg.values[[1]])
   mon.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)      
   sub.matrix <- matrix(unlist(mon.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
   rm(mon.avg.values)
   mon.avg.matrix[!flag,] <- sub.matrix
   rm(sub.matrix)
   ncvar_put(mon.ncs,varid=mon.name,vals=mon.avg.matrix,
            start=c(1,lat.ix,1),count=c(-1,1,-1))
   rm(mon.avg.matrix)
   gc()
}

##----------------------------------------------------------------------------------------------
annual.extremes.for.model <- function(ext.name,ext.ncs,lat.ix,n.lon,yearly.fac,flag,
                                      ext.list) {
      
   ##Variables
   flen <- sum(!flag)
   sub.list <- ext.list[!flag]
   rm(ext.list)
   ext.fx <- ext.fxns[[ext.name]]

   ext.avg.values <- foreach(
                             data=sub.list,
                             .export=c('yearly.fac','ext.fx')
                             ) %do% {
                                  ext.avg.values <- ext.fx(data,yearly.fac)
                             }
   rm(sub.list)                         
   ncol <- length(ext.avg.values[[1]])
   ext.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
   sub.matrix <- matrix(unlist(ext.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
   rm(ext.avg.values)
   ext.avg.matrix[!flag,] <- sub.matrix
   rm(sub.matrix)
   ncvar_put(ext.ncs,varid=ext.name,vals=ext.avg.matrix,
             start=c(1,lat.ix,1),count=c(-1,1,-1))
    rm(ext.avg.matrix)		
}


##----------------------------------------------------------------------------------------------
annual.quantiles.for.model <- function(quant.name,quant.value,quant.ncs,lat.ix,n.lon,yearly.fac,flag,
                                       quant.list) {

   flen <- sum(!flag)     
   sub.list <- quant.list[!flag] 
   quant.fx  <- function(data,fac,pctl){tapply(data,fac,quantile,pctl,na.rm=T)}

   ##Variables
   quant.val <- as.numeric(quant.value)/1000
   print(quant.val)
   quant.avg.values <- foreach(
                           data=sub.list,
                           .export=c('yearly.fac','quant.fx','quant.val')
                           ) %do% {
                                quant.avg.values <- quant.fx(data,yearly.fac,quant.val)
                           }
   ncol <- length(quant.avg.values[[1]])
   quant.avg.matrix <- matrix(NA,nrow=n.lon,ncol=ncol)
   sub.matrix <- matrix(unlist(quant.avg.values),nrow=flen,ncol=ncol,byrow=TRUE)
   rm(quant.avg.values)
   quant.avg.matrix[!flag,] <- sub.matrix
   rm(sub.matrix)
   ncvar_put(quant.ncs,varid=quant.name,vals=quant.avg.matrix,
             start=c(1,lat.ix,1),count=c(-1,1,-1))
  rm(quant.avg.matrix)		
  gc()
}

##--------------------------------------------------------------
##****************************************************************

if (1==1) {
  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }
}

##tmp.dir <- tmpdir ##'/local_temp/ssobie/prism/' ##tmpdir
tmp.dir <- paste0(tmpdir,'/',gcm,'_',type,'_',pctl,'_',varname,'/') ## '/local_temp/ssobie/prism/' ##tmpdir

##gcm <- 'CNRM-CM5'
##scenario <- 'rcp85'
##run <- 'r1i1p1'
##interval <- '1951-2100'
##type <- 'annual_quantiles'
##varname <- 'tasmin'
##pctl <- '004'

##Latitude Bands
lat.st <- seq(48.1,59.9,0.1) ##format(seq(48.0,59.9,0.1),nsmall=1)
lat.en <- seq(48.2,60.0,0.1) ##format(seq(48.1,60.0,0.1),nsmall=1)

###lat.st <- seq(48.1,48.1,0.1) ##format(seq(48.0,59.9,0.1),nsmall=1)
###lat.en <- seq(48.2,48.2,0.1) ##format(seq(48.1,60.0,0.1),nsmall=1)

len <- length(lat.st)

if (cedar) {
   base.dir <- '/scratch/ssobie/prism/'
} else {
   base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
}

data.dir <- paste0(base.dir,gcm,'/lat_split/')
##data.dir <- paste0('/storage/data/climate/downscale/CMIP5_delivery/lat_split/')

template.dir <- paste0(base.dir,gcm,'/template/',scenario,'/',type,'/')
template.file <- list.files(path=template.dir,pattern=varname)

if (grepl('annual_quantile',type)) {
 template.file <- template.file[grep(pctl,template.file)]
}


##Move data to local storage for better I/O
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}
dir.create(paste0(tmp.dir,scenario,'/',type),recursive=T)

rtm <- proc.time()
print('Copying')
print(template.file)

##file.copy(from=template.dir,to=paste0(tmp.dir,scenario,'/'),overwrite=TRUE,recursive=TRUE) ##
file.copy(from=paste0(template.dir,template.file),to=paste0(tmp.dir,scenario,'/',type,'/'),overwrite=TRUE,recursive=TRUE) ##


print('Move template time')
print(proc.time()-rtm)

##---------------------------------------------------------------------------
if (type=='annual') {
  ##Annual Average Files for writing
  print('Ann Avg opening')
  ann.name <- varname
  avg.type <- switch(ann.name,
                     pr='total',
                     tasmax='average',
                     tasmin='average')
  ann.dir <- paste0(tmp.dir,scenario,'/annual/')
  ann.file <- paste0(ann.dir,ann.name,'_annual_',avg.type,'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  ann.ncs <- nc_open(ann.file,write=TRUE)  
  common.lat <- ncvar_get(ann.ncs,'lat')
}
##---------------------------------------------------------------------------
if (type=='annual_extremes') {
##Annual Block Maxima Files for writing
  print('Ann extremes opening')
  ext.name <- varname
  ext.type <- switch(ext.name,
                     pr='maximum',tasmax='maximum',tasmin='minimum')
  ext.dir <- paste0(tmp.dir,scenario,'/annual_extremes/')
  ext.file <- paste0(ext.dir,ext.name,'_annual_',ext.type,'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  ext.ncs <- nc_open(ext.file,write=TRUE)
  common.lat <- ncvar_get(ext.ncs,'lat')
}

##---------------------------------------------------------------------------
if (grepl('annual_quantile',type)) {
  ##Annual Quantile Files for writing
  print('Ann quantiles opening')
  quant.dir <- paste0(tmp.dir,scenario,'/annual_quantiles/')
  quant.file <- paste0(quant.dir,varname,'_annual_quantile_',sprintf("%s",pctl),'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  print(quant.file)
  quant.ncs <- nc_open(quant.file,write=TRUE)
  common.lat <- ncvar_get(quant.ncs,'lat')
}

##---------------------------------------------------------------------------
if (type=='seasonal') {
##Seasonal Average Files for writing
  print('Seasonal averages opening')
  seas.name <- varname
  seas.dir <- paste0(tmp.dir,scenario,'/seasonal/')
  seas.file <- paste0(seas.dir,seas.name,'_seasonal_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  seas.ncs <- nc_open(seas.file,write=TRUE)
  common.lat <- ncvar_get(seas.ncs,'lat')
}
##---------------------------------------------------------------------------
if (type=='monthly') {
  ##Monthly Average Files for writing
  print('monthly avg opening')
  mon.name <- varname
  mon.dir <- paste0(tmp.dir,scenario,'/monthly/')
  mon.file <- paste0(mon.dir,mon.name,'_monthly_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
  mon.ncs <- nc_open(mon.file,write=TRUE)
  common.lat <- ncvar_get(mon.ncs,'lat')
}

##---------------------------------------------------------------------------
##---------------------------------------------------------------------------

##Iterate over the latitude files
for (i in 1:len) {
  lat.interval <- paste0(sprintf('%.1f',lat.st[i]),'-',sprintf('%.1f',lat.en[i]))

  data.input <- paste0(varname,'_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  print(paste0('Copy ',data.input))
  file.copy(paste0(data.dir,"/",data.input),tmp.dir,overwrite=TRUE)

  print('Data opening')
  input.file <- paste0(varname,'_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',interval,'_',lat.interval,'.nc')
  input.nc <- nc_open(paste0(tmp.dir,input.file),write=FALSE)
  input.dates <- netcdf.calendar(input.nc)
  yearly.fac <- as.factor(format(input.dates,'%Y'))
  seasonal.fac <- get.seasonal.fac(input.dates)
  monthly.fac <- as.factor(format(input.dates,'%Y-%m'))  

  lon <- ncvar_get(input.nc,'lon')
  lat <- ncvar_get(input.nc,'lat')
  n.lon <- length(lon)
  n.lat <- length(lat)
  lat.match <- which(common.lat %in% lat)
  lat.bnds <- c(lat.match[1],(tail(lat.match,1)-lat.match[1])+1)

  print('Latitude bands:')
  print(lat.bnds)

  for (j in 1:n.lat) { ##n.lon) {
    lat.ix <- lat.match[j]
    ltm <- proc.time()

    print(paste0('Latitude: ',j,' of ',n.lat))
    input.subset <- ncvar_get(input.nc,varname,start=c(1,j,1),count=c(-1,1,-1))
    flag <- is.na(input.subset[,1])
    input.list <- vector(mode='list',length=n.lon)

    rtm <- proc.time()     
    input.list <- lapply(seq_len(nrow(input.subset)), function(k) input.subset[k,])
    rm(input.subset)
    print('Lapply convert to list')
    print(proc.time()-rtm) 

    ##----------------------------------------------------------
    ##Annual Averages 
    if (type=='annual') {
    rtm <- proc.time()     
    annual.averages.for.model(ann.name,ann.ncs,lat.ix,n.lon,yearly.fac,flag,
                              input.list)
    print('Annual Averages Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    ##Seasonal Averages 
    if (type=='seasonal') {
    rtm <- proc.time()     
      seasonal.averages.for.model(seas.name,seas.ncs,lat.ix,n.lon,seasonal.fac,flag,
                                  input.list)
    print('Seasonal Averages Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    ##Monthly Averages 
    if (type=='monthly') {
    rtm <- proc.time()     
      monthly.averages.for.model(mon.name,mon.ncs,lat.ix,n.lon,monthly.fac,flag,
                                input.list)

    print('Monthly Averages Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    if (type=='annual_extremes') {
    ##Annual Extremes
    rtm <- proc.time()     
      annual.extremes.for.model(ext.name,ext.ncs,lat.ix,n.lon,yearly.fac,flag,
                                input.list)
    print('Annual Extremes Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------
    if (type=='annual_quantiles') {
    ##Annual Quantiles
    rtm <- proc.time()     
    annual.quantiles.for.model(varname,pctl,quant.ncs,lat.ix,n.lon,yearly.fac,flag,
                              input.list)
    print('Annual Quantiles Time')
    print(proc.time()-rtm) 
    }
    ##----------------------------------------------------------

    rm(input.list)

    print('Lon loop time')
    print(proc.time()-ltm)
  }##Longitude Loop
  nc_close(input.nc)

  print('Removing lat band files')
  file.remove(paste0(tmp.dir,"/",data.input))
}##Latitude File Loop

##Move back
write.dir <- paste0(base.dir,gcm)

##file.copy(from=paste0(tmp.dir,scenario,"/",type,"/"),to=paste0(write.dir,'/',scenario,'/'),overwrite=TRUE,recursive=TRUE)
print('from')
print(paste0(tmp.dir,scenario,"/",type,"/",template.file))
print('to')
print(paste0(write.dir,'/',scenario,'/',type,'/'))


  file.copy(from=paste0(tmp.dir,scenario,"/",type,"/",template.file),to=paste0(write.dir,'/',scenario,'/',type,'/'),overwrite=TRUE)

if (type=='annual') {
   nc_close(ann.ncs)
}

if (type=='seasonal') {
    nc_close(seas.ncs)
}

if (type=='monthly') {
  nc_close(mon.ncs)
}

if (type=='annual_extremes') {
    nc_close(ext.ncs)
}

if (type=='annual_quantiles') {
    nc_close(quant.ncs)
}



clean.up <- paste("rm ",tmp.dir,scenario,"/",type,"/*",sep='')
print(clean.up)
##system(clean.up)

clean.up <- paste("rmdir ",tmp.dir,scenario,"/",type,sep='')
print(clean.up)
##system(clean.up)

print("Total Elapsed Time:")
print(proc.time()-ptm)

