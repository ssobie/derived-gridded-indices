##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

space.subset <- function(input.file,output.file,bounds) {
  lons <- format(bounds$lon,nsmall=1)
  lats <- format(bounds$lat,nsmall=1)               
  work <- paste0("ncks -O -d lon,",lons[1],",",lons[2]," -d lat,",lats[1],",",lats[2]," ",input.file," ",output.file)
  print(work)
  system(work)
}

run.subset <- function(var.name,gcm,scenario,run,type,clim,subclim,bounds,
                              region,proj.dir) {

    print(gcm)
    if (grepl('ETCCDI',var.name)) {
      read.dir <- paste(proj.dir,'bccaq_gcm_bc_subset/',gcm,'/',scenario,'/climdex/climatologies/',sep='')
      write.dir  <- paste(proj.dir,'assessment_subsets/',region,'/',scenario,'/climdex/',gcm,'/',sep='')
    } else {
      read.dir <- paste(proj.dir,'bccaq_gcm_bc_subset/',gcm,'/',scenario,'/',type,'/climatologies/',sep='')
      write.dir  <- paste(proj.dir,'assessment_subsets/',region,'/',scenario,'/',type,'/',gcm,'/',sep='')
    }

    if (!file.exists(write.dir)) {
       dir.create(write.dir,recursive=T)
    }

    intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')

    for (interval in intervals) {

      input.file <- paste0(read.dir,var.name,'_',clim,'_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      output.file <- paste0(write.dir,var.name,'_',region,'_',type,'_average_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      space.subset(input.file,output.file,bounds)
    }
}

rps.subset <- function(var.name,gcm,scenario,run,rperiod,bounds,
                              region,proj.dir) {

    print(gcm)
    read.dir <- paste(proj.dir,'bccaq_gcm_bc_subset/',gcm,'/',scenario,'/return_periods/',gcm,'/',sep='')
    write.dir  <- paste(proj.dir,'assessment_subsets/',region,'/',scenario,'/return_periods/',gcm,'/',sep='')
    if (!file.exists(write.dir)) {
       dir.create(write.dir,recursive=T)
    }
    intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
    for (interval in intervals) {
      input.file <- paste0(read.dir,var.name,'_RP',rperiod,'_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      output.file <- paste0(write.dir,var.name,'_',region,'_RP',rperiod,'_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      space.subset(input.file,output.file,bounds)
    }
}


##***********************************************************************
gcms <- 'CCSM4' ##c('ACCESS1-0','CanESM2','inmcm4','CNRM-CM5','CSIRO-Mk3-6-0') ##c('ACCESS1-0','CanESM2','CNRM-CM5','inmcm4')
scenario <- 'rcp85'
run <- 'r2i1p1'

proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/'
region <- 'van_coastal_health' ##'willow_road'

##bounds <- list(lon=c(-130.0,-119.8),lat=c(53.9,60.1)) ##Northeast

bounds <- list(lon=c(-125.5,-121.5),lat=c(48.0,51.5)) ##Van Coastal Health

##bounds <- list(lon=c(-129.7,-124.8),lat=c(51.0,53.5)) ##Bella Bella Health

##bounds <- list(lon=c(-123.4,-121.0),lat=c(52.8,54.4)) ##Willow Road



##---------------------------------------------------
##Standard Indices

run.standard <- function(gcms) {
  for (gcm in gcms) {
    ##Annual
    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
                    type='annual',clim='average_annual_total',bounds=bounds,region=region,proj.dir=proj.dir) 
##    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='maximum_annual_total',bounds=bounds,region=region,proj.dir=proj.dir)  
##    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='minimum_annual_total',bounds=bounds,region=region,proj.dir=proj.dir) 
##    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='standard_deviation_annual_total',bounds=bounds,region=region,proj.dir=proj.dir) 

    run.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,
                    type='annual',clim='average_annual',bounds=bounds,region=region,proj.dir=proj.dir) 
    run.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,
                    type='annual',clim='average_annual',bounds=bounds,region=region,proj.dir=proj.dir) 

    ##Monthly
    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
                    type='monthly',clim='monthly_average',bounds=bounds,region=region,proj.dir=proj.dir) 
    run.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,
                    type='monthly',clim='monthly_average',bounds=bounds,region=region,proj.dir=proj.dir) 
    run.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,
                    type='monthly',clim='monthly_average',bounds=bounds,region=region,proj.dir=proj.dir) 

    ##Seasonal
    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
                    type='seasonal',clim='seasonal_average',bounds=bounds,region=region,proj.dir=proj.dir) 
    run.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,
                    type='seasonal',clim='seasonal_average',bounds=bounds,region=region,proj.dir=proj.dir) 
    run.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,
                    type='seasonal',clim='seasonal_average',bounds=bounds,region=region,proj.dir=proj.dir) 
  }
}

run.return.periods <- function(gcms,rperiod) {
   ##Return Periods
  for (gcm in gcms) {
##     rps.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,rperiod=rperiod,
##                  bounds=bounds,region=region,proj.dir=proj.dir) 
##     rps.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,rperiod=rperiod,
##                  bounds=bounds,region=region,proj.dir=proj.dir) 
     rps.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,rperiod=rperiod,
                  bounds=bounds,region=region,proj.dir=proj.dir) 
  }
}

run.degree.days <- function(gcms) {
   ##Degree Days
   for (gcm in gcms) {
     run.subset(var.name='cdd',gcm=gcm,scenario=scenario,run=run,
                   type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='fdd',gcm=gcm,scenario=scenario,run=run,
                  type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='gdd',gcm=gcm,scenario=scenario,run=run,
                  type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='hdd',gcm=gcm,scenario=scenario,run=run,
                  type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
  }

}


run.climdex <- function(gcms) {
   ##Climdex

##
climdex.names <- c('dtrETCCDI','rx1dayETCCDI','rx2dayETCCDI', 'rx5dayETCCDI',
                   'sdiiETCCDI','r10mmETCCDI','r20mmETCCDI','cddETCCDI','cwdETCCDI',
                   'r95pETCCDI','r99pETCCDI','prcptotETCCDI','r95daysETCCDI','r99daysETCCDI',
                    'fdETCCDI','suETCCDI','su30ETCCDI','idETCCDI','trETCCDI','txxETCCDI','txnETCCDI',
                   'tnnETCCDI','tnxETCCDI','rx2dayETCCDI','gslETCCDI')
##climdex.names <- 'gslETCCDI'
  for (gcm in gcms) {
    for (clim.name in climdex.names) {
        run.subset(var.name=clim.name,gcm=gcm,scenario=scenario,run=run,
                   type='annual',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir)    
      
    }
  }
}


