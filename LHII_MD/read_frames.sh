#! /bin/bash

#SBATCH --job-name=read_frames
#SBATCH --output=read_frames.out
#SBATCH --time=0-01:00:00
#SBATCH --ntasks-per-node=20
#SBATCH --mem=50gb

module load lang/python/anaconda/3.8-2020.07

source activate openmm

python read_frames.py

