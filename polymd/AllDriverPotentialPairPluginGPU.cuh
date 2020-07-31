// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander / Everyone is free to add additional potentials

/*! \file AllDriverPotentialPairGPU.cuh
    \brief Declares driver functions for computing all types of pair forces on the GPU
*/

#ifndef __ALL_DRIVER_POTENTIAL_PAIR_GPU_CUH__
#define __ALL_DRIVER_POTENTIAL_PAIR_GPU_CUH__

#include "hoomd/md/PotentialPairGPU.cuh"
//The above line might be potentially problematic

//! Compute lj pair forces on the GPU with PairEvaluatorLJPlugin
cudaError_t gpu_compute_ljplugintemp_forces(const pair_args_t& pair_args,
                                      const Scalar2 *d_params);
//! Compute lj pair forces on the GPU with PairEvaluatorLJPlugin
cudaError_t gpu_compute_forceshiftedljplugintemp_forces(const pair_args_t& pair_args,
                                      const Scalar2 *d_params);

//! Compute ludovic potential pair forces on the GPU with PairEvaluatorLudovic
cudaError_t gpu_compute_polydispersetemp_forces(const pair_args_t& pair_args,
                                      const Scalar3 *d_params);
//! Compute ludovic potential pair forces on the GPU with PairEvaluatorLudovic
cudaError_t gpu_compute_polydisperse_ljtemp_forces(const pair_args_t& pair_args,
                                      const Scalar3 *d_params);
cudaError_t gpu_compute_polydisperse_18temp_forces(const pair_args_t& pair_args,
                                      const Scalar3 *d_params);

cudaError_t gpu_compute_polydisperse_10temp_forces(const pair_args_t& pair_args,
                                      const Scalar3 *d_params);
#endif
