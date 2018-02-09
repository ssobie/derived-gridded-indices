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

run.subset <- function(var.name,gcm,scenario,run,type,bounds,
                              region,proj.dir) {

    print(gcm)
    read.dir <- paste(proj.dir,'bccaq_gcm_bc_subset/',gcm,'/',scenario,'/',type,'/',sep='')
    write.dir  <- paste(proj.dir,'assessment_subsets/',region,'/',scenario,'/',type,'/',gcm,'/',sep='')

    intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')

    for (interval in intervals) {

      input.file <- paste0(read.dir,var.name,'_',type,'_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      output.file <- paste0(write.dir,var.name,'_',region,'_',type,'_climatology_BCCAQ2_PRISM_',gcm,'_',scenario,'_',run,'_',interval,'.nc')
      space.subset(input.file,output.file,bounds)
    }
}


##***********************************************************************
var.name <- 'pr'
gcm <- 'CanESM2'
scenario <- 'rcp85'
run <- 'r1i1p1'
type <- 'monthly'

proj.dir <-  '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/'

region <- 'northeast'

bounds <- list(lon=c(-130.0,-119.8),lat=c(53.9,60.1))

run.subset('pr',gcm,scenario,run,type,bounds,region,proj.dir) 
run.subset('tasmax',gcm,scenario,run,type,bounds,region,proj.dir) 
run.subset('tasmin',gcm,scenario,run,type,bounds,region,proj.dir) 