#!/bin/bash
#SBATCH --job-name=1000ps_100fs_LHII
#SBATCH --output=1000ps_100fs_LHII.out
#SBATCH --error=1000ps_100fs_LHII.err
#SBATCH --time=0-5:40:00
#SBATCH --gpus=1
#SBATCH --mem=80gb

export ENV_NAME="openmm"
export INPUT_SCRIPT="${SLURM_SUBMIT_DIR}/run_openmm.py"
export INPCRD_PATH="${SLURM_SUBMIT_DIR}/LH2.inpcrd"
export PRMTOP_PATH="${SLURM_SUBMIT_DIR}/LH2.prmtop"

TMP_RESULT_DIR="${TMPDIR}/output"
RESULT_DIR="${SLURM_SUBMIT_DIR}/output"

# Job information
echo "Host:    $(hostname)"
echo "Time:    $(date)"
echo "Dir:     $(pwd)"
echo "Job ID:  ${SLURM_JOB_ID}"
echo "Nodelist:"
echo "  $(echo "${SLURM_JOB_NODELIST}" | uniq)"

module load lang/python/anaconda/3.8-2020.07

source activate ${ENV_NAME}

echo "Conda env:  ${ENV_NAME}"

echo "Python:     $(which python)"
echo "Job script: ${INPUT_SCRIPT}"
echo "Result dir: ${RESULT_DIR}"

echo -e "\n[Running job script]\n"

mkdir -p ${RESULT_DIR}

export OUTPUT_NAME="1000ps_100fs_LHII"
export LENGTH=1000     # ps
export INTEGRATOR_DT=2 # fs
export REPORTER_DT=100 # fs

time python ${INPUT_SCRIPT}

echo -e "\n[Job script finished]\n"

echo -e "Copying job script results to ${RESULT_DIR}\n"
#cp -v ${TMP_RESULT_DIR}/* ${RESULT_DIR}/
