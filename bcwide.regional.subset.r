##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

space.subset <- function(input.file,output.file,bounds) {
  lons <- format(bounds$lon,nsmall=1)
  lats <- format(bounds$lat,nsmall=1)               
  work <- paste0("ncks -O -d lon,",lons[1],",",lons[2]," -d lat,",lats[1],",",lats[2]," ",input.file," ",output.file)
  print(work)
  result <- system(work)
  if (result!=0) {
    browser()
  }
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
      ##output.file <- paste0(write.dir,var.name,'_',region,'_',type,'_average_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      output.file <- paste0(write.dir,var.name,'_',clim,'_',region,'_average_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      ##output.file <- paste0(write.dir,var.name,'_minimum_',region,'_',type,'_total_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      browser()
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



##---------------------------------------------------
##Standard Indices

run.standard <- function(gcms,runs,region,bounds) {
  for (g in seq_along(gcms)) {
    gcm <- gcms[g]
    run <- runs[g]
    ##Annual

##    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='average_annual_total',bounds=bounds,region=region,proj.dir=proj.dir) 
##    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='maximum_annual_total',bounds=bounds,region=region,proj.dir=proj.dir)  

##    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='minimum_annual_total',bounds=bounds,region=region,proj.dir=proj.dir) 
##    run.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='standard_deviation_annual_total',bounds=bounds,region=region,proj.dir=proj.dir) 

    run.subset(var.name='tas',gcm=gcm,scenario=scenario,run=run,
                    type='annual',clim='average_annual',bounds=bounds,region=region,proj.dir=proj.dir) 
##    run.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='average_annual',bounds=bounds,region=region,proj.dir=proj.dir) 
##    run.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,
##                    type='annual',clim='average_annual',bounds=bounds,region=region,proj.dir=proj.dir) 

    run.subset(var.name='tas',gcm=gcm,scenario=scenario,run=run,
                    type='monthly',clim='monthly_average',bounds=bounds,region=region,proj.dir=proj.dir) 
    run.subset(var.name='tas',gcm=gcm,scenario=scenario,run=run,
                    type='seasonal',clim='seasonal_average',bounds=bounds,region=region,proj.dir=proj.dir) 


 if (1==0) {
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
}

run.return.periods <- function(gcms,runs,rperiod,region,bounds) {
   ##Return Periods
  for (g in seq_along(gcms)) {
    gcm <- gcms[g]
    run <- runs[g]

     rps.subset(var.name='pr',gcm=gcm,scenario=scenario,run=run,rperiod=rperiod,
                  bounds=bounds,region=region,proj.dir=proj.dir) 
     rps.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,rperiod=rperiod,
                  bounds=bounds,region=region,proj.dir=proj.dir) 
     rps.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,rperiod=rperiod,
                  bounds=bounds,region=region,proj.dir=proj.dir) 
  }
}

run.degree.days <- function(gcms,runs,region,bounds) {
   ##Degree Days
  for (g in seq_along(gcms)) {
    gcm <- gcms[g]
    run <- runs[g]
     run.subset(var.name='cdd',gcm=gcm,scenario=scenario,run=run,
                   type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
##     run.subset(var.name='fdd',gcm=gcm,scenario=scenario,run=run,
##                  type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
##     run.subset(var.name='gdd',gcm=gcm,scenario=scenario,run=run,
##                  type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
##     run.subset(var.name='hdd',gcm=gcm,scenario=scenario,run=run,
##                  type='degree_days',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir) 
  }

}

   ##Climdex
run.climdex <- function(gcms,runs,region,bounds) {


climdex.names <- c('fdETCCDI','dtrETCCDI','rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI', 
                   'sdiiETCCDI','r10mmETCCDI','r20mmETCCDI',
                   'cddETCCDI','cwdETCCDI','cdd90ETCCDI','cddmaxETCCDI',##'cdd30ETCCDI',
                   'r95pETCCDI','r99pETCCDI','prcptotETCCDI','r95daysETCCDI','r99daysETCCDI',
                   'suETCCDI','su30ETCCDI','idETCCDI','trETCCDI','txxETCCDI','txnETCCDI',
                   'tnnETCCDI','tnxETCCDI','gslETCCDI')
  climdex.names <- c('su30ETCCDI','rx5dayETCCDI','r95pETCCDI','trETCCDI') 
  for (g in seq_along(gcms)) {
    gcm <- gcms[g]
    run <- runs[g]
    for (clim.name in climdex.names) {
        run.subset(var.name=clim.name,gcm=gcm,scenario=scenario,run=run,
                   type='annual',clim='annual',bounds=bounds,region=region,proj.dir=proj.dir)    
    }     
  }

  if (1==0) {
  climdex.names <- c('rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI','dtrETCCDI',    
                      'txxETCCDI','txnETCCDI','tnnETCCDI','tnxETCCDI')
  for (g in seq_along(gcms)) {
    gcm <- gcms[g]
    run <- runs[g]
    for (clim.name in climdex.names) {
          ##Monthly
          run.subset(var.name=clim.name,gcm=gcm,scenario=scenario,run=run,
                    type='monthly',clim='monthly',bounds=bounds,region=region,proj.dir=proj.dir) 
          ##Seasonal
          run.subset(var.name=clim.name,gcm=gcm,scenario=scenario,run=run,
                    type='seasonal',clim='seasonal',bounds=bounds,region=region,proj.dir=proj.dir) 
    }
  }
  }
}

run.quantiles <- function(gcms,runs,region,bounds) {
   ##Quantiles
  for (g in seq_along(gcms)) {
     gcm <- gcms[g]
     run <- runs[g]
     run.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,
                   type='annual_quantiles',clim='annual_quantile_975',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,
                   type='annual_quantiles',clim='annual_quantile_990',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='tasmax',gcm=gcm,scenario=scenario,run=run,
                   type='annual_quantiles',clim='annual_quantile_996',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,
                   type='annual_quantiles',clim='annual_quantile_004',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,
                   type='annual_quantiles',clim='annual_quantile_010',bounds=bounds,region=region,proj.dir=proj.dir) 
     run.subset(var.name='tasmin',gcm=gcm,scenario=scenario,run=run,
                   type='annual_quantiles',clim='annual_quantile_025',bounds=bounds,region=region,proj.dir=proj.dir) 

  }

}

##***********************************************************************
gcms <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4','MRI-CGCM3',
          'MIROC5','MPI-ESM-LR','CCSM4') 
runs <- c('r1i1p1','r1i1p1','r1i1p1','r1i1p1','r1i1p1','r1i1p1','r1i1p1','r1i1p1','r1i1p1',
          'r3i1p1','r3i1p1',
          'r2i1p1')
scenario <- 'rcp85'


proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/'
##region <- 'van_coastal_health' ##'willow_road'

##bounds <- list(lon=c(-130.0,-119.8),lat=c(53.9,60.1)) ##Northeast

##bounds <- list(lon=c(-125.5,-121.5),lat=c(48.0,51.5)) ##Van Coastal Health

##bounds <- list(lon=c(-129.7,-124.8),lat=c(51.0,53.5)) ##Bella Bella Health

##bounds <- list(lon=c(-123.4,-121.0),lat=c(52.8,54.4)) ##Willow Road

##'inshuck' ##
##regions <- c('bella_health','northeast','willow_road','van_coastal_health')

## bound.list <- list(list(lon=c(-122.9,-122.0),lat=c(49.6,50.5))) ##Inshuck
##bound.list <- list(list(lon=c(-129.7,-124.8),lat=c(51.0,53.5)), ##Bella Bella Health
##                   list(lon=c(-130.0,-119.8),lat=c(53.9,60.1)), ##Northeast                   
##                   list(lon=c(-123.4,-121.0),lat=c(52.8,54.4)), ##Willow Road
##                   list(lon=c(-125.5,-121.5),lat=c(48.0,51.5))) ##Van Coastal Health

##regions <- c('kootenays','central','northeast','interior_health')
##bound.list <- list(list(lon=c(-119.6,-113.9),lat=c(48.9,51.6)), ##Kootenays
##                   list(lon=c(-128.65,-118.0),lat=c(52.0,56.3)), ##Central
##                   list(lon=c(-130.0,-119.8),lat=c(53.9,60.1)), ##Northeast
##                   list(lon=c(-125.5,-113.8),lat=c(48.9,53.3))) ##Interior Health
##regions <- c('bella_health','van_coastal_health')
##bound.list <- list(list(lon=c(-129.7,-124.8),lat=c(51.0,53.5)), ##Bella Bella Health
##                   list(lon=c(-125.5,-121.5),lat=c(48.0,51.5))) ##Van Coastal Health

regions <- c('toquaht','interior_health','fraser_health')
bound.list <- list(list(lon=c(-125.88,-124.98),lat=c(48.7,49.2)), ##Toquaht
                   list(lon=c(-125.5,-113.8),lat=c(48.9,53.3)), ##Interior Health
                   list(lon=c(-123.65,-120.5),lat=c(48.7,50.35))) ##Fraser Health

##Run all together

for (i in seq_along(regions)) {
  region <- regions[i]
  bounds <- bound.list[[i]]
  run.standard(gcms,runs,region,bounds)
##  run.return.periods(gcms,runs,rperiod=5,region,bounds)
##  run.degree.days(gcms,runs,region,bounds)
##  run.climdex(gcms,runs,region,bounds)
##  run.quantiles(gcms,runs,region,bounds)
}