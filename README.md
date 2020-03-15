# **PolydisperseMD**

PolydisperseMD is a plug-in for HOOMD-Blue, a particle simulation toolkit, that implements specialized pair potentials for poly-disperse interacting particle system. The code is largely based on HOOMD-Blue's [MD's pair potentials](https://hoomd-blue.readthedocs.io/en/stable/module-md-pair.html).

The plugin is currently under *beta* stage. For now, this file will be a temporary documentation for the plugin. 

## **Contents** 

Files that come with this plugin:
 - CMakeLists.txt   : main CMake configuration file for the plugin
 - FindHOOMD.cmake  : script to find a HOOMD-Blue installation to link against
 - README.md        : This file
 - polymd           : Directory containing C++ and Python source codes that interacts with HOOMD-Blue


## **Installation Instructions**

Parts of the instructions were modified from the example plugin provided by HOOMD-Blue. See https://hoomd-blue.readthedocs.io/en/stable/developer.html for other useful information:

### **Requirements**
The requirements for installing the plugin is the same as standard HOOMD. See [HOOMD installation page](https://hoomd-blue.readthedocs.io/en/stable/installation.html) for details. 

### **Installing Plugin**

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

Here's a case example. Suppose that I'm installing the plugin in my personal workstation, where my username is 'yourusername' and the Python environment was Conda's 'base'. If all goes well, then I should see (as part of cmake's output) the following lines:
```console
-- Python output: /home/yourusername/anaconda3/envs/hoomd/lib/python3.7/site-packages/hoomd
-- Looking for a HOOMD installation at /home/yourusername/anaconda3/envs/hoomd/lib/python3.7/site-packages/hoomd
-- Found hoomd installation at /home/yourusername/anaconda3/envs/hoomd/lib/python3.7/site-packages/hoomd
-- Found HOOMD include directory: /home/yourusername/anaconda3/envs/hoomd/lib/python3.7/site-packages/hoomd/include
-- Found PythonLibs: /home/yourusername/anaconda3/envs/hoomd/lib/libpython3.7m.so
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
$ cmake ../ -DHOOMD_ROOT=/path/to/hoomd
```
where ${HOOMD_ROOT}/bin/hoomd is where the hoomd executable is installed. In the example above /path/to/hoomd is /home/yourusername/anaconda3/envs/hoomd/. 

Finally, you would compile and install the plugin:
```console
$ make -j4 install
```

---
**NOTE**

If hoomd is installed in a system directory (such as via an rpm or deb package). If you follow instructions from HOOMD's installation page, then this is highly unlikely. For completion, we will also provide instructions for this case. 

First, Delete the contents of your build folder. Set the environment variable HOOMD_PLUGINS_DIR inyour .bash_profile or .bashrc:
```console
export HOOMD_PLUGINS_DIR=${HOME}/hoomd_plugins  # as an example
```
When running cmake, you will add -DHOOMD_PLUGINS_DIR=${HOOMD_PLUGINS_DIR} to the options. Go back to your build folder now, and run:
```console
$ cmake ../ -DHOOMD_PLUGINS_DIR=${HOOMD_PLUGINS_DIR}
```

Now, 'make install' will install the plugins into ${HOOMD_PLUGINS_DIR} and hoomd, when launched, will look there
for the plugins.

---

## **Using PolydisperseMD with HOOMD's MD**

The plugin is a complement for HOOMD's MD, which means that you need to import the plugin and hoomd.md side-by-side. Everything you do will be just like running HOOMD's MD except at the part where you're defining pair potentials. For example (this script is taken from HOOMD-Blue's docs directly), 

```python
import hoomd
from hoomd import md
from hoomd import polymd # this is our plugin!
hoomd.context.initialize()

# Create a 10x10x10 simple cubic lattice of particles with type name A
hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=2.0, type_name='A'), n=10)

# Specify Lennard-Jones interactions between particle pairs
# We will use polymd for this
nl = md.nlist.cell()
lj = polymd.pair.lj_plugin(r_cut=2.5, nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)

# Integrate at constant temperature
md.integrate.mode_standard(dt=0.005)
hoomd.md.integrate.langevin(group=hoomd.group.all(), kT=1.2, seed=4)

# Run for 10,000 time steps
hoomd.run(10e3)
```

(More Instructions, coming soon . . .)

## **Developer Notes**

(More notes, coming soon . . .)
