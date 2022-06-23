#! /bin/bash
#SBATCH --job-name=ether_75000
#SBATCH --output=ether_75000.out
#SBATCH --time=0-03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=35gb

module load lang/python/anaconda/3.8-2020.07

source activate openmm

export QCORE_PATH=~/.local/src/Qcore/release/bin/qcore

export DCD_FILE=../../../Ether_MD/output/300ps_2fs_ether.dcd
export PRMTOP_FILE=../../../Ether_MD/forcefield_prep/bchla_in_diethyl_ether.prmtop

export FRAME_START=75000
export FRAME_END=80000

export OMP_NUM_THREADS=1
export MKL_THREADING_LAYER=TBB

time python ../../../src/run_ether_traj.py
