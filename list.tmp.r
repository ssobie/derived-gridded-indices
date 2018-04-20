##Script to clean up temporary directories

tmp.dir <- '/local_temp/ssobie/'
print(list.files(path=tmp.dir,recursive=T))

##file.copy(from='/local_temp/ssobie/prism/rcp85/climdex/gslETCCDI_ann_BCCAQ2-PRISM_CNRM-CM5_rcp85_r1i1p1_1951-2100.nc',
##          to='/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/CNRM-CM5/rcp85/climdex/',overwrite=T,recursive=T)


