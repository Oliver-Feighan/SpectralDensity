#!/bin/bash
#SBATCH --job-name=10ps_2fs_ether
#SBATCH --output=10ps_2fs_ether.out
#SBATCH --error=10ps_2fs_ether.err
#SBATCH --time=0-15:00:00
#SBATCH --gpus=1
#SBATCH --mem=80gb

export ENV_NAME="openmm"
export INPUT_SCRIPT="${SLURM_SUBMIT_DIR}/run_openmm.py"
export XML_PATH="${SLURM_SUBMIT_DIR}/cla_in_ether_System.xml"
export PDB_PATH="${SLURM_SUBMIT_DIR}/cla_in_ether.pdb"

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

export OUTPUT_NAME="10ps_2fs_Ether"
export LENGTH=10       # ps
export INTEGRATOR_DT=2 # fs
export REPORTER_DT=2 # fs

time python ${INPUT_SCRIPT}

echo -e "\n[Job script finished]\n"

echo -e "Copying job script results to ${RESULT_DIR}\n"
#cp -v ${TMP_RESULT_DIR}/* ${RESULT_DIR}/
