#!/bin/bash
#PBS -N robin_nicolas 
#PBS -q training
#PBS -A imf_lille-tma4280
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:15:00
#PBS -l select=2:ncpus=20:mpiprocs=8:ompthreads=4

cd $PBS_O_WORKDIR

echo $PBS_O_WORKDIR

module load gcc/6.3.0
module load openmpi/2.0.1

mpirun ./poisson 1024 4 
