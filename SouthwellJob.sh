#!/bin/bash -l

#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J Southwell
#SBATCH -C haswell
#SBATCH -A m708 

cd
. ./bash_scripts/southwell_init.sh

N=1

./DMEM_bash_scripts/SweepPar.sh ${1} $N
