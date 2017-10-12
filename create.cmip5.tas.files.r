##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)

make.tas.files <- function() {

scenarios <- c('piControl')
##'ACCESS1-0',
##'CanESM2',

gcm.list <- c('ACCESS1-0', 'bcc-csm1-1','CanESM2','CSIRO-Mk3-6-0', 'GFDL-CM3','HadGEM2-ES','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-MR',
              'ACCESS1-3',  'bcc-csm1-1-m','CCSM4','GFDL-ESM2G','HadGEM2-AO', 'inmcm4','IPSL-CM5B-LR','MIROC-ESM-CHEM','MRI-CGCM3',
              'BNU-ESM','CNRM-CM5','GFDL-ESM2M','HadGEM2-CC','IPSL-CM5A-LR','MIROC5','MPI-ESM-LR')

gcm.list <- c('ACCESS1-0','ACCESS1-3','BNU-ESM','CanESM2','CESM1-CAM5','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5B-LR',
              'MPI-ESM-LR','MPI-ESM-MR','NorESM1-M','bcc-csm1-1','bcc-csm1-1-m','CCSM4','CNRM-CM5','GFDL-CM3','inmcm4',
              'HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MRI-CGCM3')

gcm.list <- 'NorESM1-M'          

proj.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/global/control/'
for (gcm in gcm.list) {  
  for (scenario in scenarios) {
    print(gcm)
    read.dir <- paste0(proj.dir,gcm,'/')
    write.dir  <- read.dir
    
    tx.files <- list.files(path=read.dir,pattern='tasmax_day')
    tx.file <- tx.files[grep(scenario,tx.files)]

    tn.files <- list.files(path=read.dir,pattern='tasmin_day')
    tn.file <- tn.files[grep(scenario,tn.files)]
    
    tas.file <- gsub(pattern='tasmax',replacement='tas',tx.file)
    ##print(paste('cdo -s -O add ',read.dir,tx.file,' ',read.dir,tn.file,' ',write.dir,'tmp2.nc',sep=''))
    ##print(paste('cdo -s -O divc,2 ',write.dir,'tmp2.nc ',write.dir,tas.file,sep=''))
    
    system(paste('cdo -s -O add ',read.dir,tx.file,' ',read.dir,tn.file,' ',write.dir,'tmp2.nc',sep=''))
    system(paste('cdo -s -O divc,2 ',write.dir,'tmp2.nc ',write.dir,tas.file,sep=''))
    system(paste('ncrename -v tasmax,tas ',write.dir,tas.file,sep=''))
    
    system(paste('rm ',write.dir,'tmp2.nc ',sep=''))
  }
}
}


make.annual.avg.files <- function() {

scenario <- 'historical\\+rcp85'
gcm.list <- c('ACCESS1-0','ACCESS1-3','BNU-ESM','CanESM2','CSIRO-Mk3-6-0','GFDL-ESM2M','IPSL-CM5A-MR','IPSL-CM5B-LR',
              'MPI-ESM-LR','MPI-ESM-MR','NorESM1-M','bcc-csm1-1','bcc-csm1-1-m','CNRM-CM5','inmcm4',
              'HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MRI-CGCM3')

##  proj.dir <- '/storage/data/climate/downscale/BCCAQ2/CMIP5/global/control/'
proj.dir <- '/storage/data/projects/CanSISE/n2w/CMIP5/'
  
  for (gcm in gcm.list) {
    print(gcm)
    read.dir <- paste0(proj.dir,gcm,'/')
    write.dir  <- paste0(proj.dir,'annual/',gcm,'/')
    if (!file.exists(write.dir)) {
       dir.create(write.dir,recursive=TRUE)
    }
    ##tas.files <- list.files(path=read.dir,pattern='tas_day')
    tas.files <- list.files(path=read.dir,pattern='snc_LImon')
    tas.file <- tas.files[grep(scenario,tas.files)]
    write.file <- gsub('_LImon_','_ann_',tas.file)
    
    system(paste('cdo -s -O yearmean ',read.dir,tas.file,' ',write.dir,write.file,sep=''))

  }
}



