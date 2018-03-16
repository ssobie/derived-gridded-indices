
library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

var.name <- 'tasmin'
gcm <- 'ACCESS1-0'

seas.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',gcm,'/rcp85/seasonal/')
day.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',gcm,'/')

past.file <- paste0(day.dir,var.name,'_day_BCCAQ2_',gcm,'_rcp85_r1i1p1_1951-2000.nc')
past.nc <- nc_open(past.file)
past.time <- ncvar_get(past.nc,'time')
past.series <- netcdf.calendar(past.nc)

proj.file <- paste0(day.dir,var.name,'_day_BCCAQ2_',gcm,'_rcp85_r1i1p1_2001-2100.nc')
proj.nc <- nc_open(proj.file)
proj.time <- ncvar_get(proj.nc,'time')
proj.series <- netcdf.calendar(proj.nc)
full.series <- c(past.series,proj.series)

full.time <- c(past.time,proj.time)

  seas.ix <- sort(c(grep('*-01-15',full.series),
                    grep('*-04-15',full.series),
                    grep('*-07-15',full.series),
                    grep('*-10-15',full.series)))

new.time <- full.time[seas.ix]
new.series <- full.series[seas.ix]

seas.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',gcm,'/rcp85/seasonal/')
seas.file <- paste0(seas.dir,var.name,'_seasonal_BCCAQ2_PRISM_',gcm,'_rcp85_r1i1p1_1951-2100.nc')
seas.nc <- nc_open(seas.file,write=TRUE)
ncvar_put(seas.nc,'time',new.time)
nc_close(seas.nc)
nc_close(past.nc)
nc_close(proj.nc)