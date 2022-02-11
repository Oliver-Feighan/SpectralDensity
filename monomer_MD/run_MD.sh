#!/bin/bash

#SBATCH --job-name=monomer_MD
#SBATCH --output=monomer_MD.out
#SBATCH --error=monomer_MD.err

#SBATCH --time=0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

export MKL_THREADING_LAYER=TBB
export OMP_NUM_THREADS=1

~/.local/src/Qcore/release/bin/qcore -s "aimd(structure(file = 'monomer.xyz') xtb() output_steps=2 save_to_file=true output_coordinates='monomer_trajectory' n_steps=10000 thermostat(temperature=300 kelvin))"
