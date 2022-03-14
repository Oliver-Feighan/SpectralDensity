#! /bin/bash
#SBATCH --job-name=ether_10000_14999
#SBATCH --output=ether_10000_14999.out
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

module load lang/python/anaconda/3.8-2020.07

source activate openmm

export OMP_NUM_THREADS=1
export MKL_THREADING_LAYER=TBB

export QCORE_PATH=~/.local/src/Qcore/release/bin/qcore


export DCD_FILE=../../Ether_MD/output/100ps_2fs_Ether.dcd
export PDB_FILE=../../Ether_MD/cla_in_ether.pdb

export FRAME_START=10000
export FRAME_END=14999

time python ../../src/run_ether_traj.py
