##******************************************************************************
##******************************************************************************

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R')
##******************************************************************************
################################################################################

get.date.bounds <- function(nc) {
  dates <- netcdf.calendar(nc)
  yst <-  gsub('-','',format(head(dates,1),'%Y'))
  yen <-  gsub('-','',format(tail(dates,1),'%Y'))
  rv <- c(yst,yen)
  return(rv)
}

interp.bccaq <- function(var.name,gcm,scenario,interval,grid.file,years,base.dir,write.dir,tmp.dir) {

  print(paste('Interpolate Anomalies: ',gcm,', ',var.name,', ',interval,sep=''))
  gcm.dir <- paste(base.dir,gcm,sep='')  

  anoms.files <- list.files(path=gcm.dir,pattern=paste(var.name,'_anoms_BCCAQ2_',sep=''))
  scen.files <- anoms.files[grep(interval,anoms.files)]
  anoms.file <- scen.files[grep(scenario,scen.files)]
  print(anoms.file)

  ##move.to <- paste0('rsync -av ',gcm.dir,'/',anoms.file,' ',tmp.dir)
  ##print(move.to)
  ##system(move.to)

  file.copy(from=paste0(gcm.dir,'/',anoms.file),to=tmp.dir)
  print(list.files(path=tmp.dir))

  base.file <- paste0(tmp.dir,anoms.file)
  print('Before loop')
  print(list.files(path=tmp.dir))
  ##Split into 1-year files
    for (y in years) {
      print('Start of loop')
      print(list.files(path=tmp.dir))
      ytm <- proc.time()

      ##Subset to 5-Year File
      tmp.file <- gsub(pattern='_anoms_BCCAQ2_',replacement='_anoms_interp_BCCAQ2_',base.file)
      work <- paste('cdo seldate,',y,'-01-01T00:00,',y,'-12-31T23:59 ',base.file,' ',tmp.file,sep='')
      print(work)
      system(work)

      snc <- nc_open(tmp.file)
      start.time <- ncvar_get(snc,'time')
      dates <- get.date.bounds(snc)
      nc_close(snc)

      ##interp.file <- gsub(pattern='[0-9]{4}-[0-9]{4}',replacement=paste0(dates[1],'-',dates[2]),tmp.file)
      interp.file <- gsub(pattern='[0-9]{4}-[0-9]{4}',replacement=y,tmp.file)
      work <- paste('cdo -s remapbil,',grid.file,' ',tmp.file,' ',interp.file,sep='')
      print(work)
      system(work)

      ##move.back <- paste0('rsync -av ',interp.file,' ',gcm.dir,'/interpolated/')
      ##print(move.back)
      ##system(move.back)
      file.copy(from=interp.file,to=write.dir)

      ##clean.up <- paste0('rm ',interp.file)
      ##print(clean.up)
      ##system(clean.up)	
      file.remove(interp.file)                          
      
      ##clean.up <- paste0('rm ',tmp.file)
      ##print(clean.up)
      ##system(clean.up)	
      file.remove(tmp.file)
      print('End of loop')
      print(list.files(path=tmp.dir))

    }
    file.remove(base.file)  
    ##clean.up <- paste0('rm ',base.file) ##This is the anomaly file
    ##print(clean.up)
    ##system(clean.up)	

  gc()  
}


################################################################################
##******************************************************************************

run.interp <- function() {

  ##Requires 1981-2010 (or equivalent base period) in the /baseline directory to create anomalies from the full period
  ## (usually 1950-2000 and 2001-2100). Both of these are created using extract.bccaq.gcm.r
  ##Also need the PRISM climatologies (also using extract.bccaq.gcm.r).

  args <- commandArgs(trailingOnly=TRUE)
  for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
  }

##gcm <- 'HadGEM2-ES'
##varname <- 'pr'
##scenario <- 'rcp85'
##tmpdir <- '/local_temp/ssobie/interpolate/'

  tmp.dir <- paste0(tmpdir,gcm,'/',varname,'/')
  if (!file.exists(tmp.dir)) {
     dir.create(tmp.dir,recursive=T)
  }

##  base.dir <- '/scratch/ssobie/prism/'
##  grid.file <- '/home/ssobie/assessments/bc.prism.grid.txt'

  base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
  write.dir <- paste0('/storage/data/climate/downscale/CMIP5_delivery/',gcm,'/interpolated/')
  grid.file <- '/storage/home/ssobie/grid_files/bc.prism.grid.txt'

  years <- seq(1951,2000,by=1)	
  interp.bccaq(varname,gcm,scenario,'1951-2000',grid.file,years,base.dir,write.dir,tmp.dir)
  years <- seq(2001,2100,by=1)	
  interp.bccaq(varname,gcm,scenario,'2001-2100',grid.file,years,base.dir,write.dir,tmp.dir) 

  clean.up <- paste("rm ",tmp.dir,"/*nc" ,sep='')
  print(clean.up)
  system(clean.up)

}

run.interp()



