#!/bin/bash                                                                                                                                         
gcm="MIROC5"
run='r3i1p1'
scenario='rcp85'

qsub -N "${gcm}.an.pr" -v gcm=$gcm,run=$run,scenario=$scenario,varname="pr" run.anomalies.pbs
qsub -N "${gcm}.an.tx" -v gcm=$gcm,run=$run,scenario=$scenario,varname="tasmax" run.anomalies.pbs
qsub -N "${gcm}.an.tn" -v gcm=$gcm,run=$run,scenario=$scenario,varname="tasmin" run.anomalies.pbs

