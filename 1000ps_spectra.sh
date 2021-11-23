#!/bin/bash
#SBATCH --job-name=1000ps_analysis
#SBATCH --output=1000ps_analysis.out
#SBATCH --time=0-00:10:00
#SBATCH --mem=50gb

module load lang/python/anaconda/3.8-2020.07

source activate openmm

python analysis.py --inpcrd LHII_MD/LH2.inpcrd --prmtop LHII_MD/LH2.prmtop --dcd LHII_MD/output/1000ps_2fs_LHII.dcd --sample_period 100
