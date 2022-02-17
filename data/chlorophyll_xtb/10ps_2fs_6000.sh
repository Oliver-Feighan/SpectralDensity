#! /bin/bash

#SBATCH --job-name=10ps_2fs
#SBATCH --output=10ps_2fs_6000.out
#SBATCH --time=0-08:00:00
#SBATCH --ntasks-per-node=20
#SBATCH --mem=50G

module load lang/python/anaconda/3.8-2020.07

source activate openmm

which python

export QCORE_PATH=~/.local/src/Qcore/release/bin/qcore
export DCD_FILE=../../LHII_MD/output/10ps_2fs_LHII.dcd
export PRMTOP_FILE=../../LHII_MD/LH2.prmtop

export FRAME_START=6000
export FRAME_END=7000

export OMP_NUM_THREADS=1
export MKL_THREADING_LAYER=TBB

python ../../src/chl_xtb.py 

