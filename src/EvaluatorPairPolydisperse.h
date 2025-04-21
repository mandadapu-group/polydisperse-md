// Original code:
// Copyright (c) 2009-2024 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Modification:
// Copyright (c) 2025, Sanggeun Song, University of California, Berkeley.

#ifndef __PAIR_EVALUATOR_POLYDISPERSE_H__
#define __PAIR_EVALUATOR_POLYDISPERSE_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPolydisperse.h
    \brief Defines the pair evaluator class for the polydisperse potential
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
#ifdef __HIPCC__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif

namespace hoomd
    {
namespace md
    {

class EvaluatorPairPolydisperse
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar v0;
        Scalar eps;
        Scalar scaledr_cut;

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : v0(1.0), eps(0), scaledr_cut(2.5) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            v0 = v["v0"].cast<Scalar>();
            eps = v["eps"].cast<Scalar>();
            scaledr_cut = v["scaledr_cut"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["v0"] = v0;
            v["eps"] = eps;
            v["scaledr_cut"] = scaledr_cut;
            return v;
            }
#endif
        }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairPolydisperse(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), v0(_params.v0), eps(_params.eps), scaledr_cut(_params.scaledr_cut)
        {
        c0 = Scalar(-28.0) * v0 / pow(scaledr_cut, 12);
        c1 = Scalar(48.0) * v0 / pow(scaledr_cut, 14);
        c2 = Scalar(-21.0) * v0 / pow(scaledr_cut, 16);
        }

    // IMPORTANT NOTE: Since the radius is not provided v4, use the charge data as a substitute for diameter.
    DEVICE static bool needsCharge()
        {
        return true;
        }
    //! Accept the optional charge value
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) 
        {
        d_i = qi;
        d_j = qj;
        }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param pair_eng Output parameter to write the computed pair energy
        \param energy_shift If true, the potential must be shifted so that
        V(r) is continuous at the cutoff
        \note There is no need to check if rsq < rcutsq in this method.
        Cutoff tests are performed in PotentialPair.

        \return True if they are evaluated or false if they are not because
        we are beyond the cutoff
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        Scalar sigma = 0.5 * (d_i + d_j) * (1 - eps * fabs(d_i - d_j));
        Scalar actualcutsq = scaledr_cut * scaledr_cut * sigma * sigma;
        // compute the force divided by r in force_divr
        if (rsq < actualcutsq && v0 != 0)
            {
            Scalar r2inv = sigma * sigma * Scalar(1.0) / rsq;
            Scalar _rsq = Scalar(1.0) * rsq / (sigma * sigma);
            Scalar r6inv = r2inv * r2inv * r2inv;
            force_divr = (Scalar(12.0) * v0 * r2inv * r6inv * r6inv - Scalar(2.0) * c1
                          - Scalar(4.0) * c2 * _rsq)
                         / (sigma * sigma);

            // No energy shift is needed
            pair_eng = v0 * r6inv * r6inv + c0 + c1 * _rsq + c2 * _rsq * _rsq;

            return true;
            }
        else
            return false;
        }

    //! Polydisperse doesn't eval LRC integrals
    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

    //! Polydisperse doesn't eval LRC integrals
    DEVICE Scalar evalEnergyLRCIntegral()
        {
        return 0;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("polydisperse_pair");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar d_i;
    Scalar d_j;

    Scalar v0;
    Scalar sigma;
    Scalar eps;
    Scalar scaledr_cut;

    // Additional parameters to be compute
    Scalar c0;
    Scalar c1;
    Scalar c2;
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_POLYDISPERSE_H__
