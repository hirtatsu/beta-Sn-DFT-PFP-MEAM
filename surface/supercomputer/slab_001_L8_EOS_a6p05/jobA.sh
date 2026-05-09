#!/bin/bash
#------- qsub option -----------
#PBS -P NIFS26KISM025
#PBS -q A_S
#PBS -l select=1:ncpus=256:mpiprocs=128:ompthreads=2:mem=752gb
#PBS -l walltime=12:00:00

#------- Program execution -----------
cd ${PBS_O_WORKDIR}

# これ以降の全ての標準出力と標準エラーを log.txt に出力する
exec > log.txt 2>&1

module load intel
module load intelmpi

echo ""
echo "################################"
echo "##      Job Information       ##"
echo "################################"
echo "Job ID: ${PBS_JOBID}"
echo "Job Name: ${PBS_JOBNAME}"
echo "Queue: ${PBS_QUEUE}"
echo "Work Directory: ${PBS_O_WORKDIR}"
echo ""

mpirun -np 128 ../openmx slab_001_L8_a6p05.dat2 -nt 2

echo ""
echo "################################"
echo "##        Job Finished        ##"
echo "################################"
