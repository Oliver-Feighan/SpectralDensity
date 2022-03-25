#! /bin/bash

#SBATCH --job-name=5000ps_10fs
#SBATCH --output=5000ps_10fs_468000.out
#SBATCH --time=0-20:00:00
#SBATCH --ntasks-per-node=20
#SBATCH --mem=50G

module load lang/python/anaconda/3.8-2020.07

source activate openmm

which python

export QCORE_PATH=~/.local/src/Qcore/release/bin/qcore
export DCD_FILE=../../../LHII_MD/output/chl_5000ps_10fs_LHII.dcd
export PRMTOP_FILE=../../../LHII_MD/LH2.prmtop
export ONLY_CHL=True

export FRAME_START=468000
export FRAME_END=472000

export OMP_NUM_THREADS=1
export MKL_THREADING_LAYER=TBB

python ../../../src/chl_xtb.py 

