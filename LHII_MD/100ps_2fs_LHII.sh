#!/bin/bash
#SBATCH --job-name=openmm
#SBATCH --output=100ps_2fs_LHII.out
#SBATCH --error=100ps_2fs_LHII.err
#SBATCH --time=0-15:00:00
#SBATCH --gpus=1
#SBATCH --mem=80gb

export ENV_NAME="openmm"
export INPUT_SCRIPT="${SLURM_SUBMIT_DIR}/run_openmm_simulation.py"
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

export OUTPUT_NAME="100ps_2fs_LHII"
export LENGTH=100
export INTEGRATOR_DT=2
export REPORTER_DT=10

time python - <<EOF
import os
from pathlib import Path 

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

import time
from datetime import timedelta

print(f"Running MD with GPU platform\n")

# Files
output_name = os.environ["OUTPUT_NAME"]
results_dir_path = Path.cwd() / "output"
print(f"output folder made at {results_dir_path}")
results_dir_path.mkdir(exist_ok=True)

prmtop_path = Path(os.environ["PRMTOP_PATH"])
inpcrd_path = Path(os.environ["INPCRD_PATH"])
assert prmtop_path.exists()
assert inpcrd_path.exists()
print(f"Reading system coordinates, topology and forcefield from {inpcrd_path}, {prmtop_path}")
prmtop = AmberPrmtopFile(prmtop_path)
inpcrd = AmberInpcrdFile(inpcrd_path)

# System Configuration
print("setup system")
nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
ewaldErrorTolerance = 0.0005
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001

# Integration Options
print("setup integration")
dt = int(os.environ["INTEGRATOR_DT"]) * femtoseconds
temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 25

# Simulation Options
print("setup simulation")
steps = 1e3 * int(os.environ["LENGTH"])
equilibrationSteps = 30000
platform = Platform.getPlatformByName('OpenCL')
platformProperties = {'Precision': 'double'}
dcdReporter = DCDReporter(f'{results_dir_path}/{output_name}.dcd', int(float(os.environ["REPORTER_DT"])))
dataReporter = StateDataReporter(f'{results_dir_path}/{output_name}.txt', int(os.environ["REPORTER_DT"]), totalSteps=steps,
	step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter(f'{results_dir_path}/{output_name}.chk', 5000)

# Prepare the Simulation

print('Building system...')
topology = prmtop.topology
positions = inpcrd.positions
system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, 
	constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platformProperties=platformProperties)
	
simulation.context.setPositions(positions)
if inpcrd.boxVectors is not None:
	simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
	
# Minimize and Equilibrate

print('Performing energy minimization...')
minimize_energy_start = time.time()
simulation.minimizeEnergy()
print(f"time for energy minimization: {str(timedelta(seconds = time.time() - minimize_energy_start))}")
print()
print(f"Equilibrating: {equilibrationSteps} steps at {dt} fs time steps")
equilibration_start = time.time()
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)
print(f"time for equilibration: {str(timedelta(seconds = time.time() - equilibration_start))}")
print()

# Simulate

simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
	
print(f"Production: {steps} steps at {dt} fs time steps - total of {os.environ['LENGTH']} ps")
print(f" saving step every {int(os.environ['REPORTER_DT']) * dt} fs - total of {1e3 * int(os.environ['LENGTH'])/(dt * int(os.environ['REPORTER_DT']))} frames") 
production_start = time.time()
simulation.currentStep = 0
simulation.step(steps)
print(f"time for production: {str(timedelta(seconds = time.time() - production_start))}")
EOF

echo -e "\n[Job script finished]\n"

echo -e "Copying job script results to ${RESULT_DIR}\n"
#cp -v ${TMP_RESULT_DIR}/* ${RESULT_DIR}/
