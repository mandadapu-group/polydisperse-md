# **PolydisperseMD**

PolydisperseMD is a plug-in for HOOMD-Blue, a particle simulation toolkit, that implements specialized pair potentials for poly-disperse interacting particle system. The code is largely based on HOOMD-Blue's [MD's pair potentials](https://hoomd-blue.readthedocs.io/en/stable/module-md-pair.html).

The plugin is currently under *beta* stage.  

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


## **How to Use PolydisperseMD**

The plugin is a complement for HOOMD's MD, which means that you need to import the plugin and hoomd.md side-by-side. 

(More Instructions, coming soon . . .)

## **Developer Notes**

(More notes, coming soon . . .)
