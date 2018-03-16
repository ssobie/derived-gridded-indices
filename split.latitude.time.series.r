##******************************************************************************
##******************************************************************************

library(ncdf4)
library(PCICt)

##******************************************************************************
################################################################################

split.and.concatenate <- function(var.name,gcm,scenario,run,
				  lat.st,lat.en,year.list,base.dir,tmp.dir) {

  ##Split off one latitude band from each of the year files				  					    
  split.files <- rep('A',length=length(year.list))
  for (y in seq_along(year.list)) {
    year <- year.list[y]
    
    print(paste('Daily PRISM: ',gcm,', ',var.name,', ',year,sep=''))
    gcm.dir <- paste(base.dir,gcm,'/',sep='')
    bccaq.file <- paste0(var.name,'_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',year,'.nc')
    ##rtm <- proc.time()
    ##file.copy(from=paste0(gcm.dir,'daily_prism/',bccaq.file), to=tmp.dir)
    bccaq.tmp <- paste0(gcm.dir,'daily_prism/',bccaq.file)
          
    split.tmp <- paste0(gcm.dir,'daily_prism/',var.name,'_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',year,'_',lat.st,'-',lat.en,'.nc')
    split.files[y] <- split.tmp
    print('Split File')		
    print(split.tmp)
    stm <- proc.time()
    cmd <- paste0('ncks -O -d lat,',lat.st,',',lat.en,' ',bccaq.tmp,' ',split.tmp,sep='')
    print(cmd)
    system(cmd)
    print('Split time:')
    print(proc.time()-stm)    		
    ##file.remove(paste0('rm ',bccaq.tmp))
  }
  yst <- head(year.list,1)
  yen <- tail(year.list,1)
  
  series.file <- paste0(var.name,'_gcm_prism_BCCAQ2_',gcm,'_',scenario,'_',run,'_',yst,'-',yen,'_',lat.st,'-',lat.en,'.nc')
  ##Concatenate the Split files into one latitude band with a full time series
  atm <- proc.time()
  cat.files <- paste0('ncrcat ',paste(split.files,collapse=' '),' ',gcm.dir,'lat_split/',series.file) 
  print(cat.files)
  system(cat.files)
  print('Concatenate time:')
  print(proc.time()-atm)    		

  ##mtm <- proc.time()
  ##file.copy(from=paste0(tmp.dir,series.file),to=paste0(gcm.dir,'lat_split/'))
  file.remove(paste0(gcm.dir,'daily_prism/',split.files))
  ##file.remove(paste0(tmp.dir,'/',series.file))
  gc()

}

################################################################################
##******************************************************************************

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

##gcm <- 'CSIRO-Mk3-6-0'
##scenario <- 'rcp85'
##tmpdir <- '/local_scratch/ssobie/split/'
##varname <- 'tasmax'

tmp.dir <- tmpdir
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=T)
}

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

var.list <- varname ##c('tasmax') ##,'tasmin','pr')

lat.st <- format(seq(48.0,59.9,0.1),nsmall=1)
lat.en <- format(seq(48.1,60.0,0.1),nsmall=1)

lat.st <- format(seq(59.2,59.9,0.1),nsmall=1) ##59.9,0.1),nsmall=1)
lat.en <- format(seq(59.3,60.0,0.1),nsmall=1) ##60.0,0.1),nsmall=1)
len <- length(lat.st)

for (var.name in var.list) {
  print(var.name)

  ##-----------------------------------------
  ##Find the run and the year bounds
  check.files <- list.files(path=paste0(base.dir,gcm),
		            pattern=paste0(var.name,'_day_BCCAQ2_',gcm,'_',scenario,'*'),full.name=T)
  check.file <- check.files[grep('2001-2100',check.files)]
  print(check.file)
  file.split <- strsplit(check.file,'_')[[1]]
  run <- file.split[grep('r*i1p1',file.split)]
  print(run)

  nc <- nc_open(check.file)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  past.values <- ncvar_get(nc,'time')
  full.series <- format(past.origin + past.values*86400,'%Y')
  tst <- 1951 ##as.numeric(head(full.series,1))
  ten <- as.numeric(tail(full.series,1))
  print(paste0('Year Interval: ',tst,'-',ten))
  year.list <- seq(tst,ten,1)
  nc_close(nc)

  ##-----------------------------------------

  year.list <- seq(tst,ten,by=1)
  for (i in 1:len) {
    itm <- proc.time()

    print(lat.st[i])
    print(lat.en[i])
    split.and.concatenate(var.name,gcm,scenario,run,
                          lat.st[i],lat.en[i],year.list,base.dir,tmp.dir)
    print('Splitting loop time for one lat band:')
    print(proc.time()-itm)    		
			  

  }
}  

##clean.up <- paste("rm ",tmp.dir,"/*nc" ,sep='')
##print(clean.up)
##system(clean.up)





