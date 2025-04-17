// Original code:
// Copyright (c) 2009-2024 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modifications:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

// Include the defined classes that are to be exported to python
#include "EvaluatorPairPolydisperse.h"

#include "hoomd/md/PotentialPair.h"
#include <pybind11/pybind11.h>
#ifdef ENABLE_HIP
#include "hoomd/md/PotentialPairGPU.h"
#endif

namespace hoomd
    {
namespace md
    {

// specify the python module. Note that the name must explicitly match the PROJECT() name provided
// in CMakeLists (with an underscore in front)
PYBIND11_MODULE(_polymd, m)
    {
    detail::export_PotentialPair<EvaluatorPairPolydisperse>(m, "PotentialPairPolydisperse");
#ifdef ENABLE_HIP
    detail::export_PotentialPairGPU<EvaluatorPairPolydisperse>(m, "PotentialPairPolydisperseGPU");
#endif
    }

    } // end namespace md
    } // end namespace hoomd
