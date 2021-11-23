#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --output=test_openmm.out
#SBATCH --error=test_openmm.err
#SBATCH --time=0-08:00:00
#SBATCH --gpus=1
#SBATCH --mem=30gb

module load lang/python/anaconda/3.8-2020.07

ENV_NAME="openmm"

source activate ${ENV_NAME}

# Job information
echo "Host:    $(hostname)"
echo "Time:    $(date)"
echo "Dir:     $(pwd)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "Nodelist:"
echo "  $(echo "${SLURM_JOB_NODELIST}" | uniq)"

echo "Conda env:  ${ENV_NAME}"
echo "Python:     $(which python)"

python -m openmm.testInstallation
python -m simtk.testInstallation
