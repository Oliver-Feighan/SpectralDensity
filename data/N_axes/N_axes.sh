#! /bin/bash

#SBATCH --job-name=N_axes_300ps_2fs
#SBATCH --output=N_axes_300ps_2fs_0.out
#SBATCH --time=0-0:05:00
#SBATCH --ntasks-per-node=5
#SBATCH --mem=30G

module load lang/python/anaconda/3.8-2020.07

source activate openmm

which python

export DCD_FILE=../../LHII_MD/output/chl_300ps_2fs_LHII.dcd
export TOP_FILE=../../LHII_MD/LH2.prmtop

time python ../../src/N_axes.py 
