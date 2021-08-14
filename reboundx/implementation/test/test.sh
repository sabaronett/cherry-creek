#!/bin/sh
#PBS -j oe
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:01:00
#PBS -m abe
#PBS -M barons2@unlv.nevada.edu
#PBS -N test_$1
#PBS -q small

cd $PBS_O_WORKDIR
module load conda
conda activate rebx-3.4.1
echo Test no. $1
