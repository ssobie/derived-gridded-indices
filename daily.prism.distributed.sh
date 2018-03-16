#!/bin/bash                                                                                                                                         
gcm='CNRM-CM5'
scenario='rcp85'
varname='tasmin'

qsub -v gcm=$gcm,varname=$varname,scenario=$scenario,tmp='first' run.daily.prism.pbs
qsub -v gcm=$gcm,varname=$varname,scenario=$scenario,tmp='second' run.daily.prism.pbs
qsub -v gcm=$gcm,varname=$varname,scenario=$scenario,tmp='third' run.daily.prism.pbs
qsub -v gcm=$gcm,varname=$varname,scenario=$scenario,tmp='fourth' run.daily.prism.pbs
qsub -v gcm=$gcm,varname=$varname,scenario=$scenario,tmp='fifth' run.daily.prism.pbs


