#!/bin/bash

#SBATCH --job-name=hessian
#SBATCH --output=hessian_run.out
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G

~/.local/src/Qcore/release/bin/qcore -f json hessian.in > hessian.out
