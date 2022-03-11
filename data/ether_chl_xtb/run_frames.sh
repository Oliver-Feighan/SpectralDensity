#! /bin/bash
#SBATCH --job-name=ether_0_4999
#SBATCH --output=ether_0_4999
#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

#module load lang/python/anaconda/3.8-2020.07

#source activate openmm

export QCORE_PATH=~/qcore/cmake-build-release/bin/qcore

export DCD_FILE=../../Ether_MD/output/10ps_2fs_Ether.dcd
export PDB_FILE=../../Ether_MD/cla_in_ether.pdb

export FRAME_START=0
export FRAME_END=5

time python ../../src/run_ether_traj.py
