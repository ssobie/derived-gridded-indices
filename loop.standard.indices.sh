#!/bin/bash                                                                                                                                         
gcm="CNRM-CM5"
run='r1i1p1'
scenario='rcp85'

##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="pr" run.separate.derived.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="tasmax" run.separate.derived.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='annual',varname="tasmin" run.separate.derived.pbs

##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='seasonal',varname="pr" run.separate.derived.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='seasonal',varname="tasmax" run.separate.derived.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='seasonal',varname="tasmin" run.separate.derived.pbs

##qsub -N "${gcm}.mon.pr" -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="pr" run.separate.derived.pbs
##qsub -N "${gcm}.mon.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="tasmax" run.separate.derived.pbs
##qsub -N "${gcm}.mon.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='monthly',varname="tasmin" run.separate.derived.pbs

##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_extremes',varname="pr" run.separate.derived.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_extremes',varname="tasmax" run.separate.derived.pbs
##qsub -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_extremes',varname="tasmin" run.separate.derived.pbs


qsub -N "${gcm}.qt.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmax",pctl="975" run.separate.derived.pbs
qsub -N "${gcm}.qt.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmax",pctl="990" run.separate.derived.pbs
qsub -N "${gcm}.qt.tx" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmax",pctl="996" run.separate.derived.pbs

##qsub -N "${gcm}.qt.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmin",pctl="004" run.separate.derived.pbs
##qsub -N "${gcm}.qt.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmin",pctl="010" run.separate.derived.pbs
qsub -N "${gcm}.qt.tn" -v gcm=$gcm,run=$run,scenario=$scenario,type='annual_quantiles',varname="tasmin",pctl="025" run.separate.derived.pbs
