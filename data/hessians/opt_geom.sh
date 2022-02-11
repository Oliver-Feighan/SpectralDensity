#!/bin/bash
#SBATCH --job-name=opt_geom
#SBATCH --output=opt_geom_run.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=0-01:00:00
#SBATCH --mem=5G


export MKL_THREADING_LAYER=TBB
export OMP_NUM_THREADS=20

~/.local/src/Qcore/release/bin/qcore opt_geom.in > opt_geom.out
