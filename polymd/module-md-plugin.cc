// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander All developers are free to add the calls needed to export their modules
#include "AllPluginPairPotentials.h"
#include "hoomd/md/PotentialPair.h"

// include GPU classes
#ifdef ENABLE_CUDA
#include "hoomd/md/PotentialPairGPU.h"
#endif

#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
namespace py = pybind11;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
    create the hoomd python module and define the exports here.
*/
PYBIND11_MODULE(_polymd, m)
    {
    export_PotentialPair<PotentialPairLJPlugin>(m, "PotentialPairLJPlugin");
    export_PotentialPair<PotentialPairForceShiftedLJPlugin>(m, "PotentialPairForceShiftedLJPlugin");
    export_PotentialPair<PotentialPairPolydisperse>(m, "PotentialPairPolydisperse");
    export_PotentialPair<PotentialPairPolydisperseLJ>(m, "PotentialPairPolydisperseLJ");
    export_PotentialPair<PotentialPairPolydisperse18>(m, "PotentialPairPolydisperse18");
    export_PotentialPair<PotentialPairPolydisperse10>(m, "PotentialPairPolydisperse10");
    export_PotentialPair<PotentialPairPolydisperseLJ106>(m, "PotentialPairPolydisperseLJ106");

#ifdef ENABLE_CUDA
    export_PotentialPairGPU<PotentialPairLJPluginGPU, PotentialPairLJPlugin>(m, "PotentialPairLJPluginGPU");
    export_PotentialPairGPU<PotentialPairForceShiftedLJPluginGPU, PotentialPairForceShiftedLJPlugin>(m, "PotentialPairForceShiftedLJPluginGPU");
    export_PotentialPairGPU<PotentialPairPolydisperseGPU, PotentialPairPolydisperse>(m, "PotentialPairPolydisperseGPU");
    export_PotentialPairGPU<PotentialPairPolydisperseLJGPU, PotentialPairPolydisperseLJ>(m, "PotentialPairPolydisperseLJGPU");
    export_PotentialPairGPU<PotentialPairPolydisperse18GPU, PotentialPairPolydisperse18>(m, "PotentialPairPolydisperse18GPU");
    export_PotentialPairGPU<PotentialPairPolydisperse10GPU, PotentialPairPolydisperse10>(m, "PotentialPairPolydisperse10GPU");
    export_PotentialPairGPU<PotentialPairPolydisperseLJ106GPU, PotentialPairPolydisperseLJ106>(m, "PotentialPairPolydisperseLJ106GPU");
#endif
    }
