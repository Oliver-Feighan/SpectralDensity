#! /bin/bash

#SBATCH --job-name=100ps_10fs
#SBATCH --output=100ps_10fs.out
#SBATCH --time=0-01:00:00
#SBATCH --nodes=20
#SBATCH --mem=10G

module load lang/python/anaconda/3.8-2020.07

source activate openmm

which python

export QCORE_PATH=~/.local/src/Qcore/release/bin/qcore
export DCD_FILE=../LHII_MD/output/1ps_2fs.dcd
export PRMTOP_FILE=../LHII_MD/LH2.prmtop

python chl_xtb.py 

