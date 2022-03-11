#! /bin/bash
#SBATCH --job-name=ether_0_4999
#SBATCH --output=ether_0_4999
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

module load lang/python/anaconda/3.8-2020.07

source activate openmm
