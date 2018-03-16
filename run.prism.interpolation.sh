#!/bin/bash 
#SBATCH --time=01-00:00
#SBATCH --nodes=1
#SBATCH --mem=20000M
#SBATCH --job-name=accprm45tx
#SBATCH --output=acc_tx_45.out 

module load r/3.3.3
module load cdo

gcm='ACCESS1-0'
scenario='rcp45'
varname='tasmax'

cd /home/ssobie/assessments/
echo "Current working directory is `pwd`"

nohup R CMD BATCH "--args tmpdir='$SLURM_TMPDIR' gcm='$gcm' scenario='$scenario' varname='$varname'" interpolate.anomalies.r "./${gcm}.${scenario}.interp.${varname}.out"
