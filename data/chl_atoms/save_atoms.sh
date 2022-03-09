#! /bin/bash

#SBATCH --job-name=10ps_2fs
#SBATCH --output=10ps_2fs_0.out
#SBATCH --time=0-00:10:00
#SBATCH --ntasks-per-node=20
#SBATCH --mem=50G

module load lang/python/anaconda/3.8-2020.07

source activate openmm

which python

python ../../src/save_chl.py

