#!/bin/bash
#------- qsub option -----------
#PBS -P NIFS26KISM025
#PBS -q A_S
#PBS -l select=1:ncpus=256:mpiprocs=128:ompthreads=2:mem=752gb
#PBS -l walltime=24:00:00

#------- Program execution -----------
cd ${PBS_O_WORKDIR}

exec > log.txt 2>&1

module load intel
module load intelmpi

echo "Job ID: ${PBS_JOBID}"

mpirun -np 128 ../openmx slab_001_L12_Relax_recon.dat -nt 2

echo "Job Finished"
