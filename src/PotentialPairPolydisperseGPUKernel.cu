// Original code:
// Copyright (c) 2009-2024 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025 Sanggeun Song, University of California, Berkeley.

#include "EvaluatorPairPolydisperse.h"
#include "hoomd/md/PotentialPairGPU.cuh"

namespace hoomd
    {
namespace md
    {
namespace kernel
    {
template __attribute__((visibility("default"))) hipError_t
gpu_compute_pair_forces<EvaluatorPairPolydisperse>(const pair_args_t& pair_args,
                                              const EvaluatorPairPolydisperse::param_type* d_params);
    } // end namespace kernel
    } // end namespace md
    } // end namespace hoomd
