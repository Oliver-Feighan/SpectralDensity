#!/bin/bash

#SBATCH --job-name=run_frames_7500_9999
#SBATCH --output=run_frames_7500_9999.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --time=0-03:00:00
#SBATCH --mem=15gb


export OMP_NUM_THREADS=1
export MKL_THREADING_LAYER=TBB

hostname

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "Job started"
echo "$dt"

module load lang/python/anaconda/3.8-2020.07

source activate openmm

export QCORE_PATH=~/.local/src/Qcore/release/bin/qcore
export FRAME_START=7500
export FRAME_END=9999
export TRAJ_FILE=truncated_trajectory.xyz 
export N_ATOMS=79

time python ../../src/run_xyz_traj.py
