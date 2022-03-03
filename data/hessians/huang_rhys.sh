#!/bin/bash

#SBATCH --job-name=Huang-Rhys
#SBATCH --output=Huang-Rhys.out
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=0-03:00:00

hostname

module load lang/python/anaconda/3.8-2020.07

source activate openmm

export OMP_NUM_THREADS=1
export MKL_THREADING_LAYER=TBB

export MOLDEN_FILE=opt_truncated.molden
export XYZ_FILE=opt_truncated.xyz

export QCORE_PATH=~/.local/src/Qcore/release/bin/qcore

time python ../../src/huang_rhys.py
