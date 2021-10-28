# Author: Matthew Hancock
# Date: 10/28/2021
# Description: Run an MD simulation using the OpenMM MD engine and AMBER force field to generate decoy structures from the ubiquitin native structure. The simulated tempering algorithm is used to explore the conformational space. The simulation output is saved to an output DCD file along with log data from the simulation.
from pathlib import Path
import time
from openmm.app import *
from openmm import *
from openmm.unit import *


interval = 1000
data_dir = Path(Path.home(), "xtal_benchmark/decoys/data")
tmp_dir = Path(Path.home(), "xtal_benchmark/decoys/tmp")
prmtop = AmberPrmtopFile(str(Path(data_dir, "input/topology.top")))
inpcrd = AmberInpcrdFile(str(Path(data_dir, "input/coords.crd")))
traj_file = Path(data_dir, "output/1ubq_st.dcd")
log_file = Path(data_dir, "output/1ubq_log.txt")
st_log_file = Path(data_dir, "output/1ubq_st_log.txt")

# Build the system from AMBERTools files
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# Create reporters.
dcd_reporter = DCDReporter(
    file=str(traj_file),
    reportInterval=interval
)

state_reporter = StateDataReporter(
    file=str(log_file),
    reportInterval=interval,
    step=True,
    potentialEnergy=True,
    kineticEnergy=True,
    totalEnergy=True,
    temperature=True,
    volume=True,
    density=True,
    speed=True
)

# Create main + 2 equilibration integrators.
integrator_100 = LangevinMiddleIntegrator(
    100*kelvin,
    1/picosecond,
    0.002*picoseconds
)
integrator_200 = LangevinMiddleIntegrator(
    200*kelvin,
    1/picosecond,
    0.002*picoseconds
)
integrator = LangevinMiddleIntegrator(
    300*kelvin,
    1/picosecond,
    0.002*picoseconds
)

# Equilibrate simulation at 100K.
simulation_100 = Simulation(
    prmtop.topology,
    system,
    integrator_100
)
simulation_100.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation_100.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
simulation_100.reporters.append(state_reporter)
print("100K equilibration")
simulation_100.minimizeEnergy()
simulation_100.step(10000)
simulation_100.saveState(str(Path(tmp_dir, "sim_100.state")))

# Equilibrate simulation at 200K.
simulation_200 = Simulation(
    prmtop.topology,
    system,
    integrator_200
)
simulation_200.reporters.append(state_reporter)
print("200K equilibration")
simulation_200.loadState(str(Path(tmp_dir, "sim_100.state")))
simulation_200.minimizeEnergy()
simulation_200.step(10000)
simulation_200.saveState(str(Path(tmp_dir, "sim_200.state")))

# Perform simulated tempering between 300-450K.
simulation = Simulation(
    prmtop.topology,
    system,
    integrator
)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)
print("Begin simulated annealing")
simulation.loadState(str(Path(tmp_dir, "sim_200.state")))
simulation.minimizeEnergy()

simulated_tempering = SimulatedTempering(
    simulation,
    numTemperatures=100,
    minTemperature=300*kelvin,
    maxTemperature=450*kelvin,
    reportInterval=interval,
    reportFile=str(st_log_file)
)

# Run for a microsecond (10E-6)
t1 = time.time()
simulated_tempering.step(50000000)
# simulated_tempering.step(50000)
t2 = time.time()
print("TIME: {}".format(t2-t1))

# simulation.step(200000)
