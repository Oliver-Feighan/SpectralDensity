#!/bin/bash
#SBATCH --job-name=hessian
#SBATCH --output=hessian_run.out
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=5G

export MKL_THREADING_LAYER=TBB
export OMP_NUM_THREADS=1

~/.local/src/Qcore/release/bin/qcore -s "hess := hessian(structure(file='${MONOMER}') xtb() save_normal_modes = '${MONOMER/xyz/molden}')" > ${MONOMER/xyz/out}
