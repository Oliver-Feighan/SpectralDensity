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


# Integration Options
print("setup integration")
dt = int(os.environ["INTEGRATOR_DT"])
dt_fs = dt * femtoseconds

print(f"Integrator timestep : {dt_fs}")

reporter_dt = int(os.environ["REPORTER_DT"]) # (int) time between saving frames in the reporters
reporter_dt_fs = reporter_dt * femtosecond   # (openmm.unit.femtosecond)  "
reporter_dt_timesteps = reporter_dt / dt     # time between saving frames in units of timesteps - saving every 10fs at 2fs dt should give saving every 5 timesteps



print(f"Reporter timestep : {reporter_dt_fs}")


nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
ewaldErrorTolerance = 0.0005
constraints = HBonds
rigidWater = True


topology = prmtop.topology
positions = inpcrd.positions

system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, 
                             constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)

temperature = 300*kelvin
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 25
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

integrator = LangevinMiddleIntegrator(temperature, friction, dt)

constraintTolerance = 0.000001
integrator.setConstraintTolerance(constraintTolerance)

platform = Platform.getPlatformByName('OpenCL')
platformProperties = {'Precision': 'double'}
simulation = Simulation(topology, system, integrator, platformProperties=platformProperties)

simulation.context.setPositions(positions)
if inpcrd.boxVectors is not None:
simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# Minimize and Equilibrate
print('Performing energy minimization...')
minimize_energy_start = time.time()
simulation.minimizeEnergy()
print(f"time for energy minimization: {str(timedelta(seconds = time.time() - minimize_energy_start))}\n")


equilibrationSteps = 30000
print(f"Equilibrating: {equilibrationSteps} steps at {dt} fs time steps")

equilibration_start = time.time()

simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

print(f"Time for equilibration: {str(timedelta(seconds = time.time() - equilibration_start))}\n")


# Simulate
steps = 1e3 * int(os.environ["LENGTH"]) / dt

dcdReporter = DCDReporter(f'{results_dir_path}/{output_name}.dcd', reporter_dt_timesteps)

dataReporter = StateDataReporter(f'{results_dir_path}/{output_name}.txt', reporter_dt_timesteps, totalSteps=steps,
                                 step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')

checkpointReporter = CheckpointReporter(f'{results_dir_path}/{output_name}.chk', 5000)

simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)

print(f"Production: {steps} steps at {dt} fs time steps - total of {(steps * dt_fs).in_units_of(picosecond)} ps")
print(f" saving step every {reporter_dt_fs} - total of {steps / reporter_dt_timesteps} frames") 
production_start = time.time()
simulation.currentStep = 0
simulation.step(steps)
print(f"time for production: {str(timedelta(seconds = time.time() - production_start))}")