#!/bin/bash                                                                                                                                         
gcm="inmcm4"
run='r1i1p1'
scenario='rcp85'
rperiod='20'

##qsub -N "${gcm}.rp${rperiod}.pr" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="pr",interval='1971-2000' run.bcwide.return.periods.pbs
##qsub -N "${gcm}.rp${rperiod}.pr" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="pr",interval='2011-2040' run.bcwide.return.periods.pbs
##qsub -N "${gcm}.rp${rperiod}.pr" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="pr",interval='2041-2070' run.bcwide.return.periods.pbs
##qsub -N "${gcm}.rp${rperiod}.pr" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="pr",interval='2071-2100' run.bcwide.return.periods.pbs

qsub -N "${gcm}.rp${rperiod}.tasmax" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmax",interval='1971-2000' run.bcwide.return.periods.pbs
qsub -N "${gcm}.rp${rperiod}.tasmax" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmax",interval='2011-2040' run.bcwide.return.periods.pbs
qsub -N "${gcm}.rp${rperiod}.tasmax" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmax",interval='2041-2070' run.bcwide.return.periods.pbs
qsub -N "${gcm}.rp${rperiod}.tasmax" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmax",interval='2071-2100' run.bcwide.return.periods.pbs

qsub -N "${gcm}.rp${rperiod}.tasmin" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmin",interval='1971-2000' run.bcwide.return.periods.pbs
qsub -N "${gcm}.rp${rperiod}.tasmin" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmin",interval='2011-2040' run.bcwide.return.periods.pbs
qsub -N "${gcm}.rp${rperiod}.tasmin" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmin",interval='2041-2070' run.bcwide.return.periods.pbs
qsub -N "${gcm}.rp${rperiod}.tasmin" -v gcm=$gcm,run=$run,scenario=$scenario,rperiod=$rperiod,varname="tasmin",interval='2071-2100' run.bcwide.return.periods.pbs



