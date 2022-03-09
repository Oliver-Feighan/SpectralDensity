#!/bin/bash

#SBATCH --job-name=truncated_MD
#SBATCH --output=truncated_MD.out
#SBATCH --error=truncated_MD.err

#SBATCH --time=0-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

export MKL_THREADING_LAYER=TBB
export OMP_NUM_THREADS=1

~/.local/src/Qcore/release/bin/qcore -s "aimd(structure(file = '../hessians/opt_truncated.xyz') xtb() output_steps=2 save_to_file=true output_coordinates='truncated_trajectory' output_energies='truncated_energies' n_steps=100000 thermostat(temperature=300 kelvin))"
