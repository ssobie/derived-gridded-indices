#!/bin/bash                                                                                                                                         
gcm="CNRM-CM5"
run='r1i1p1'
scenario='rcp85'

##qsub -N "clim.su" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="su" run.separate.climdex.pbs
##qsub -N "clim.su30" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="su30" run.separate.climdex.pbs
##qsub -N "clim.id" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="id" run.separate.climdex.pbs
##qsub -N "clim.fd" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="fd" run.separate.climdex.pbs
##qsub -N "clim.tr" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="tr" run.separate.climdex.pbs

##qsub -N "clim.r10mm" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="r10mm" run.separate.climdex.pbs
##qsub -N "clim.r20mm" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="r20mm" run.separate.climdex.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="sdii" run.separate.climdex.pbs
##qsub -N "clim.prcp" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="prcptot" run.separate.climdex.pbs
qsub -N "clim.cwd" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="cwd" run.separate.climdex.pbs


##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="rx1day" run.separate.climdex.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="rx2day" run.separate.climdex.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="rx5day" run.separate.climdex.pbs

##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="txx" run.separate.climdex.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="txn" run.separate.climdex.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="tnn" run.separate.climdex.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="tnx" run.separate.climdex.pbs

