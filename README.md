
# **PolydisperseMD**

PolydisperseMD is a plug-in for HOOMD-Blue, a particle simulation toolkit, that implements specialized pair potentials for poly-disperse interacting particle system. The code is largely based on HOOMD-Blue's [MD's pair potentials](https://hoomd-blue.readthedocs.io/en/stable/module-md-pair.html).

The plugin is ready to use! For now, this file will be a temporary documentation for the plugin. 

## **Contents** 

Files that come with this plugin:
 - CMakeLists.txt   : main CMake configuration file for the plugin
 - FindHOOMD.cmake  : script to find a HOOMD-Blue installation to link against
 - README.md        : This file
 - polymd           : Directory containing C++ and Python source codes that interacts with HOOMD-Blue


## **Installation Instructions**

Parts of the instructions were modified from the example plugin provided by HOOMD-Blue. See https://hoomd-blue.readthedocs.io/en/stable/developer.html for other useful information:

### Step 1: **Check Requirements**
The requirements for installing the plugin is the same as standard HOOMD. See [HOOMD installation page](https://hoomd-blue.readthedocs.io/en/stable/installation.html) for details. 

### Step 2: **Install  Plugin**

The process is similar to installing HOOMD.  First, git clone the project:
```console
$ git clone https://github.com/mandadapu-group/polydisperse-md
```

Next, configure your build.
```console
$ cd polydisperse-md
$ mkdir build
$ cd build
$ cmake ../ -DENABLE_MPI=ON -DENABLE_CUDA=ON
```

In this step, CMake will try to find the usual required packages (including LLVM). However, it will also try to find a HOOMD installation. Check your CMake output! 

Here's a case example. Suppose that I'm installing the plugin in my personal workstation, where my username is `yourusername` and the Python environment was a virtual environment named `iluvbase`. If all goes well, then I should see (as part of CMake's output) the following lines:
```console
-- Python output: /home/yourusername/anaconda3/envs/iluvbase/lib/python3.7/site-packages/hoomd
-- Looking for a HOOMD installation at /home/yourusername/anaconda3/envs/iluvbase/lib/python3.7/site-packages/hoomd
-- Found hoomd installation at /home/yourusername/anaconda3/envs/iluvbase/lib/python3.7/site-packages/hoomd
-- Found HOOMD include directory: /home/yourusername/anaconda3/envs/iluvbase/lib/python3.7/site-packages/hoomd/include
-- Found PythonLibs: /home/yourusername/anaconda3/envs/iluvbase/lib/libpython3.7m.so
```

If not, then the following output could be found:
```console
CMake Error at FindHOOMD.cmake:46 (message):
  Could not find hoomd installation, either set HOOMD_ROOT or set
  PYTHON_EXECUTABLE to a python which can find hoomd
Call Stack (most recent call first):
  CMakeLists.txt:8 (include)
```
if the above message is what you found, then delete the contents of your build folder. Next, re-run cmake with the following build option:
```
$ cmake ../ -DENABLE_MPI=ON -DHOOMD_ROOT=/path/to/hoomd
```
where `${HOOMD_ROOT}/bin/hoomd` is where the hoomd executable is installed. In the example above `/path/to/hoomd` is `/home/yourusername/anaconda3/envs/hoomd/`. 

Finally, you would compile and install the plugin:
```console
$ make -j4 install
```
## **Additional Notes**

If hoomd is installed in a system directory (such as via an rpm or deb package). If you follow instructions from HOOMD's installation page, then this is highly unlikely. For completion, we will also provide instructions for this case. 

First, Delete the contents of your build folder. Set the environment variable `HOOMD_PLUGINS_DIR` in `your .bash_profile`or  `.bashrc`:
```console
export HOOMD_PLUGINS_DIR=${HOME}/hoomd_plugins  # as an example
```
When running CMake, you will add `-DHOOMD_PLUGINS_DIR=${HOOMD_PLUGINS_DIR}` to the options. Go back to your build folder now, and run:
```console
$ cmake ../ -DENABLE_MPI=ON -DHOOMD_PLUGINS_DIR=${HOOMD_PLUGINS_DIR}
```

Now, `make install` will install the plugins into `${HOOMD_PLUGINS_DIR}`. When hoomd is launched, it will look in that directory for the plugins. 

---

## **Using PolydisperseMD with HOOMD's MD**

The plugin is a complement for HOOMD's MD, which means that you need to import the plugin and hoomd.md side-by-side. Everything you do will be just like running HOOMD's MD except at the part where you're defining pair potentials. For example, (this script is taken from HOOMD-Blue's docs directly), 

```python
import hoomd
import hoomd.md as md
import hoomd.polymd as polymd
import numpy as np
from numpy.random import uniform, seed

seed(0)
hoomd.context.initialize("--mode=cpu --notice-level=2");

#Set up "equilibration" and "production" runs
deltat = 0.002
totalsteps = 10/deltat
eqtotalsteps = 10/deltat
snapshots = 1000
kT = 0.25

#Initialize Our Own Configuration using a Snapshot
rho = 1.00
LParticles = 64;
NParticles = LParticles**2
dmax = 1.0
dmin = 0.5

Length = LParticles
MyBox = hoomd.data.boxdim(L=Length, dimensions=2)
snap = hoomd.data.make_snapshot(N=NParticles, box=MyBox, particle_types=['A'])
snap.particles.types = ['A']

def placePolydisperseOnSquare(snap):
    for i in range(LParticles):
        for j in range(LParticles):
            snap.particles.position[i*LParticles+j,0] = Length*(i/LParticles-0.5)
            snap.particles.position[i*LParticles+j,1] = Length*(j/LParticles-0.5)
            snap.particles.position[i*LParticles+j,2] = 0
            snap.particles.diameter[i*LParticles+j] = uniform(dmin,dmax)

placePolydisperseOnSquare(snap)
system = hoomd.init.read_snapshot(snap);
nl = md.nlist.cell()

#Set up the pair potential
poly12 = polymd.pair.polydisperse(r_cut=4.0,nlist=nl,model='polydisperse12')
poly12.pair_coeff.set('A', 'A',v0=1.0,eps=0.2,scaledr_cut=1.25)
poly12.set_params(mode="no_shift")


#Set up NVT thermostat
all = hoomd.group.all();
md.integrate.mode_standard(dt=deltat);
integrator = md.integrate.nvt(group=all, kT=kT, tau=50*deltat);
integrator.randomize_velocities(seed=339021)

#Equilibration Run on NVT Ensemble
eqsamplingtime = int(eqtotalsteps/100)
if eqsamplingtime == 0:
    eqsamplingtime = 1
logger = hoomd.analyze.log(filename='mylogeq1.log', period=eqsamplingtime, quantities=['temperature','potential_energy','kinetic_energy','momentum'], phase=0)
hoomd.run_upto(eqtotalsteps);

#Production Run on NVE Ensemble
integrator.disable()
integrator = md.integrate.nve(group=all)

samplingtime = int(totalsteps/snapshots)
if samplingtime == 0:
    samplingtime = 1

hoomd.dump.gsd(filename="restart1.gsd", group=hoomd.group.all(), dynamic=['attribute','momentum'],truncate=True, period=samplingtime, phase=0)
hoomd.dump.gsd(filename="dump1.gsd", group=hoomd.group.all(), dynamic=['attribute','momentum'],period=samplingtime, phase=0)
logger = hoomd.analyze.log(filename='mylog1.log', period=samplingtime, quantities=['temperature','pressure_xx','pressure_yy','pressure_xy','potential_energy','kinetic_energy','momentum'],phase=0)
hoomd.run_upto(eqtotalsteps+totalsteps);
```

Note that the pair potentials available in this plugin are limited to a particular class where the pair potential is given by:

![equation](https://latex.codecogs.com/gif.latex?%5Cphi%28r/%5Csigma_%7B%5Calpha%20%5Cbeta%7D%29%20%3D%20%5Cbegin%7Bcases%7D%20v_0%20%5Cleft%5B%5Cleft%28%5Cdfrac%7B%5Csigma_%7B%5Calpha%20%5Cbeta%7D%7D%7Br%7D%5Cright%29%5Em-%5Cleft%28%5Cdfrac%7B%5Csigma_%7B%5Calpha%20%5Cbeta%7D%7D%7Br%7D%5Cright%29%5En%5Cright%5D&plus;%5Csum_%7Bk%3D0%7D%5Eq%20c_k%20%5Cleft%28%5Cfrac%7Br%5E%7B%5Calpha%20%5Cbeta%7D%7D%7B%5Csigma_%7B%5Calpha%20%5Cbeta%7D%7D%20%5Cright%20%29%5E%7B2k%7D%26%20r/%5Csigma_%7B%5Calpha%20%5Cbeta%7D%20%5Cleq%20%5Ctilde%7Br%7D_c%20%5C%5C%200%20%26%20%5Ctext%7Botherwise%7D%20%5Cend%7Bcases%7D)

![equation](https://latex.codecogs.com/gif.latex?%5Csigma_%7B%5Calpha%20%5Cbeta%7D%20%3D%20%5Cfrac%7B1%7D%7B2%7D%5Cleft%28%5Csigma_%5Calpha%20&plus;%5Csigma_%5Cbeta%5Cright%29%281-%5Cvarepsilon%7C%5Csigma_%5Calpha%20-%20%5Csigma_%5Cbeta%7C%29)

The first term in the first equation is the standard repulsive and attractive interaction. The second term is an even polynomial ensuring smoothness up to q-th order at the cut off radius. 

The models available to use and currently implemented are as follows:

|   Model Name      |   q       |   m       |   n       | 
|   :--------       |   :--:    |   :--:    |   :--:    |
|   polydisperse12  |   2       |   12      |   0       |
|   polydisperse18  |   2       |   18      |   0       |
|   polydisperselj  |   2       |   12      |   6       |
|   polydisperse10  |   3       |   10      |   0       |
|   polydisperse106 |   2       |   10      |   6       |

You will see in polymd/pair.py file that there are other pair potentials, but I haven't thoroughly tested them or haven't checked their implementation in a long time! So be please be aware. 

(More Instructions, coming soon . . .)

## **Developer Notes**

(More notes, coming soon . . .)
