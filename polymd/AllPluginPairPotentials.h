// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander / Anyone is free to add their own pair potentials here

#ifndef __PAIR_POTENTIALS_PLUGIN__H__
#define __PAIR_POTENTIALS_PLUGIN__H__

#include "hoomd/md/PotentialPair.h"
#include "EvaluatorPairLJPlugin.h"
#include "EvaluatorPairForceShiftedLJPlugin.h"
#include "EvaluatorPairPolydisperse.h"
#include "EvaluatorPairPolydisperseLJ.h"
#include "EvaluatorPairPolydisperse18.h"
#include "EvaluatorPairPolydisperse10.h"
#include "EvaluatorPairPolydisperseLJ106.h"
#include "EvaluatorPairPolydisperseYukawa.h"

#ifdef ENABLE_CUDA
#include "hoomd/md/PotentialPairGPU.h"
#include "AllDriverPotentialPairPluginGPU.cuh"
#endif

/*! \file AllPairPotentials.h
    \brief Handy list of typedefs for all of the templated pair potentials in hoomd
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

//! Pair potential force compute for lj forces
typedef PotentialPair<EvaluatorPairLJPlugin> PotentialPairLJPlugin;
typedef PotentialPair<EvaluatorPairForceShiftedLJPlugin> PotentialPairForceShiftedLJPlugin;
typedef PotentialPair<EvaluatorPairPolydisperse> PotentialPairPolydisperse;
typedef PotentialPair<EvaluatorPairPolydisperseLJ> PotentialPairPolydisperseLJ;
typedef PotentialPair<EvaluatorPairPolydisperse18> PotentialPairPolydisperse18;
typedef PotentialPair<EvaluatorPairPolydisperse10> PotentialPairPolydisperse10;
typedef PotentialPair<EvaluatorPairPolydisperseLJ106> PotentialPairPolydisperseLJ106;
typedef PotentialPair<EvaluatorPairPolydisperseYukawa> PotentialPairPolydisperseYukawa;

#ifdef ENABLE_CUDA
//! Pair potential force compute for lj forces on the GPU
typedef PotentialPairGPU< EvaluatorPairLJPlugin, gpu_compute_ljplugintemp_forces > PotentialPairLJPluginGPU;
typedef PotentialPairGPU< EvaluatorPairForceShiftedLJPlugin, gpu_compute_forceshiftedljplugintemp_forces > PotentialPairForceShiftedLJPluginGPU;
typedef PotentialPairGPU< EvaluatorPairPolydisperse, gpu_compute_polydispersetemp_forces > PotentialPairPolydisperseGPU;
typedef PotentialPairGPU< EvaluatorPairPolydisperseLJ, gpu_compute_polydisperse_ljtemp_forces > PotentialPairPolydisperseLJGPU;
typedef PotentialPairGPU< EvaluatorPairPolydisperse18, gpu_compute_polydisperse_18temp_forces > PotentialPairPolydisperse18GPU;
typedef PotentialPairGPU< EvaluatorPairPolydisperse10, gpu_compute_polydisperse_10temp_forces > PotentialPairPolydisperse10GPU;
typedef PotentialPairGPU< EvaluatorPairPolydisperseLJ106, gpu_compute_polydisperse_LJ106temp_forces > PotentialPairPolydisperseLJ106GPU;
typedef PotentialPairGPU< EvaluatorPairPolydisperseYukawa, gpu_compute_polydisperse_yukawatemp_forces > PotentialPairPolydisperseYukawaGPU;
#endif

#endif // __PAIR_POTENTIALS_PLUGIN_H__
