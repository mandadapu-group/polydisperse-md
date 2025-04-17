
# **PolydisperseMD**

PolydisperseMD is a plug-in for HOOMD-Blue, a particle simulation toolkit, that implements specialized pair potentials for poly-disperse interacting particle system. The code is largely based on HOOMD-Blue's [MD's pair potentials](https://hoomd-blue.readthedocs.io/en/stable/hoomd/md/module-pair.html).

The plugin is ready to use and is compatible with HOOMD-Blue v4 and v5.

## **Contents** 

Files that come with this plugin:
 - CMakeLists.txt   : main CMake configuration file for the plugin
 - README.md        : This file
 - src              : Directory containing C++ and Python source codes that interacts with HOOMD-Blue


## **Installation Instructions**

Parts of the instructions were modified from the hoomd-component-template repository provided by HOOMD-Blue. See https://hoomd-blue.readthedocs.io/en/stable/components.html for other useful information.

### Step 1: **Check Requirements**
The requirements for installing the plugin is the same as standard HOOMD. See [HOOMD installation page, building from source](https://hoomd-blue.readthedocs.io/en/stable/building.html) for details. 

### Step 2: **Install  Plugin**

The process is similar to installing HOOMD.  First, git clone the project:
```console
$ git clone https://github.com/mandadapu-group/polydisperse-md
```

Next, configure your build and install.
```console
$ cd polydisperse-md
$ cmake -B build -S .
$ cmake --build ./build
$ cmake --install ./build
```

In this step, CMake will try to find the usual required packages. However, it will also try to find a HOOMD installation. Check your CMake output! 

## **Using PolydisperseMD with HOOMD's MD**

The plugin is a complement for HOOMD's MD, which means that you need to import the plugin and hoomd.md side-by-side. Everything you do will be just like running HOOMD's MD except at the part where you're defining pair potentials. For example,  

```python
import hoomd
import hoomd.md as md
import hoomd.polymd as polymd
import numpy as np
from numpy.random import uniform
import gsd.hoomd

seed(0)
simulation = hoomd.Simulation(device = hoomd.device.CPU(notice_level = 2), seed = 0)

# Set up the parameter for molecular dynamics
deltat = 0.002
totalsteps = 10/deltat
snapshots = 1000
kT = 0.25

# Set up the filename
temp_gsd_filename = 'MDrun_temp.gsd'
result_gsd_filename = 'MDrun.gsd'
log_filename = 'MDrun.log'

# Initialize Our Own Configuration using a Snapshot
rho = 1.00
LParticles = 64
NParticles = LParticles**2
dmax = 1.0
dmin = 0.5

def placePolydisperseOnSquare(position, diameter):
    for i in range(LParticles):
        for j in range(LParticles):
            position[i*LParticles+j,0] = Length*(i/LParticles-0.5)
            position[i*LParticles+j,1] = Length*(j/LParticles-0.5)
            position[i*LParticles+j,2] = 0
            diameter[i*LParticles+j] = uniform(dmin,dmax)

frame = gsd.hoomd.Frame()
frame.particles.N = NParticles
position = np.zeros((NParticles, 3))
diameter = np.zeros((NParticles,))
placePolydisperseOnSquare(position, diameter)
frame.particles.position = position
frame.particles.typeid = [0] * NParticles
frame.particles.types = ['A']
frame.configuration.box = [LParticles, LParticles, 0, 0, 0, 0]
# Because HOOMD v4 or higher does not provide diameter in evaluating potential, we use charge instead of diameter!
frame.particles.charge = diameter

# Define integrator and neighboring list
integrator = md.Integrator(dt = deltat)
nl = md.nlist.Cell(buffer = 0.4)

# Set up the pair potential
poly12 = polymd.pair.PolydispersePair(nlist = nl, default_r_cut = 4.0)
poly12.params[('A', 'A')] = dict(v0 = 1.0, eps = 0.2, scaledr_cut = 1.25)
integrator.forces.append(poly12)

# Compute thermodynamic properties
thermodynamic_properties = md.compute.ThermodynamicQuantities(filter = hoomd.filter.All())
simulation.operations.computes.append(thermodynamic_properties)

# Set the Logger and gsd writer
samplingtime = int(totalsteps / snapshots)
if samplingtime == 0:
    samplingtime = 1
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(simulation, quantities=['timestep'])
logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 'kinetic_energy', 'potential_energy'])
trigger = hoomd.trigger.Periodic(period = samplingtime, phase = 0)
gsd_writer = hoomd.write.GSD(trigger = trigger, filename = temp_gsd_filename, filter = hoomd.filter.All(), mode='wb', 
                            dynamic=['particles/position', 'particles/charge', 'particles/velocity'])
simulation.operations.writers.append(gsd_writer)

# Show the status of simulation
class Status:
    def __init__(self, simulation):
        self.simulation = simulation

    @property
    def seconds_remaining(self):
        try:
            return (
                self.simulation.final_timestep - self.simulation.timestep
            ) / self.simulation.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining)).split('.',2)[0]

log_print = 1000
status_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
status_logger.add(simulation, quantities=['timestep', 'tps'])
status = Status(simulation)
status_logger[('Status', 'etr')] = (status, 'etr', 'string')
table_status = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=log_print), logger=status_logger)
simulation.operations.writers.append(table_status)

# Set up NVT thermostat
with open(log_filename, 'w') as f:
    table = hoomd.write.Table(trigger = trigger, logger = logger, output = f)
    simulation.operations.writers.append(table)
    simulation.run(totalsteps, write_at_start=True)
    gsd_writer.flush()

# Convert charge to diameter in GSD file
data = gsd.hoomd.open(temp_gsd_filename, 'r')
NParticles = data[0].particles.N
with gsd.hoomd.open(result_gsd_filename, 'w') as f:
    nframe = len(data)
    newframe = gsd.hoomd.Frame()
    for i in range(0, nframe):
        frame = data[i]
        newframe = frame
        newframe.particles.diameter = frame.particles.charge
        newframe.particles.charge = np.zeros((NParticles,))
        f.append(newframe)
os.remove(temp_gsd_filename)

```

Note that the pair potential in this plugin is given by:

```math
\phi(r_{ij} / \sigma_{ij}) = v_0  \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{12} + c_0 + c_2 \left( \frac{r_{ij}}{\sigma_{ij}} \right)^2 + c_4 \left( \frac{r_{ij}}{\sigma_{ij}} \right)^4
```
```math
$\sigma_{ij} = 0.5(\sigma_i + \sigma_j)(1 - \epsilon | \sigma_i - \sigma_j|)$,
```
where the constants $c_0 = -28 v_0 / \tilde{r}_c^{12}$, $c_2 = 48 v_0 / \tilde{r}_c^{14}$, and $c_4 = -21 v_0 / \tilde{r}_c^{16}$.

## **Developer Notes**

To make your own potential, you may need to modify `EvaluatorPairPolydisperse.h` file located in `src` directory. Specifically, you should uopdate the `EvaluatorPairPolydisperse` and `evalForceAndEnergy` functions as needed. If additional variables are required, you can define them as `Scalar` type at the bottom of the header file.
