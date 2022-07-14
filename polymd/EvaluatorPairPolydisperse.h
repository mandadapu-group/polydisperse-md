// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_POLYDISPERSE_H__
#define __PAIR_EVALUATOR_POLYDISPERSE_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairPolydisperse.h
    \brief Defines the pair evaluator class for Polydisperse potentials
    \details As the prototypical example of a MD pair potential, this also serves as the primary
   documentation and base reference for the implementation of pair evaluators.
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
    //! Class for evaluating the Polydisperse pair potential
    /*! <b>General Overview</b>
    
     TO DO: Write a general overview

    */
        class EvaluatorPairPolydisperse
        {
            public:
                //! Define the parameter type used by this pair potential evaluator
                struct param_type
                {
                    Scalar v0;
                    Scalar rcut;
                    Scalar eps;
                    Scalar c0;
                    Scalar c1;
                    Scalar c2;
                    Scalar m_expnt;
                    Scalar n_expnt;

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
                    param_type() : v0(0.0), m_expnt(0), n_expnt(0), rcut(0.0), eps(0.0), c0(0.0), c1(0.0), c2(0.0)
                    {}

                    param_type(pybind11::dict v, bool managed = false)
                    {
                        v0 = v["v0"].cast<Scalar>();
                        rcut = v["rcut"].cast<Scalar>();
                        eps = v["eps"] .cast<Scalar>();
                        c0 = v["c0"].cast<Scalar>();
                        c1 = v["c1"].cast<Scalar>();
                        c2 = v["c2"].cast<Scalar>();
                        m_expnt = v["m_expnt"].cast<Scalar>();
                        n_expnt = v["n_expnt"].cast<Scalar>();
                    }

                    pybind11::dict asDict()
                    {
                        pybind11::dict v;
                        v["eps"] = eps;//#0.0;
                        v["v0"] = v0;//#0.0;
                        v["m_expnt"] = m_expnt;//0;
                        v["n_expnt"] = n_expnt;
                        v["rcut"] = rcut;
                        v["c2"] = c2;
                        v["c1"] = c1;
                        v["c0"] = c0;
                        return v;
                    }
#endif
                }
#ifdef SINGLE_PRECISION
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
                : rsq(_rsq), rcutsq(_rcutsq), v0(_params.v0), eps(_params.eps), rcut(_params.rcut), m_expnt(static_cast<int>(_params.m_expnt)), n_expnt(static_cast<int>(_params.n_expnt)), mhalf(static_cast<int>(_params.m_expnt/2)), nhalf(static_cast<int>(_params.n_expnt/2)), c0(_params.c0), c1(_params.c1), c2(_params.c2)
                {
                }

            //! Polydisperse particles obviously use diameter
            DEVICE static bool needsDiameter()
                {
                return true;
                }
            //! Accept the optional diameter values
            /*! \param di Diameter of particle i
                \param dj Diameter of particle j
            */
            DEVICE void setDiameter(Scalar di, Scalar dj) 
                {
                   d_i = di;
                   d_j = dj; 
                }

            //! Polydisperse doesn't use charge
            DEVICE static bool needsCharge()
                {
                return false;
                }
            //! Accept the optional diameter values
            /*! \param qi Charge of particle i
                \param qj Charge of particle j
            */
            DEVICE void setCharge(Scalar qi, Scalar qj) { }

            //! Evaluate the force and energy
            DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
            {
                Scalar sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
                Scalar sigmasq = sigma*sigma;
                Scalar actualcutsq = rcut*rcut*sigmasq;

                // compute the force divided by r in force_divr
                if (rsq <= actualcutsq && v0 != 0)
                {
                    Scalar r2inv = sigmasq*Scalar(1.0)/rsq;
                    Scalar r2 = Scalar(1.0)*rsq/(sigmasq);
                    Scalar rattr = 1.0;
                    Scalar rrep = 1.0;
                    
                    //Figure out how many times to for loop to get the attractive and repulsive part of the potentials
                    for (int i = 0; i < mhalf; i++)
                    {
                        rrep *= r2;
                        if (i < nhalf)
                        {
                            rattr *= r2;
                        }
                    }
                    force_divr = (Scalar(m_expnt)*v0*r2inv*rrep-Scalar(n_expnt)*v0*r2inv*rattr-Scalar(2.0)*c1 -Scalar(4.0)*c2*rsq)/(sigmasq);
                    
                    //No energy shift is needed
                    pair_eng = v0*(rrep-rattr)+c0+c1*rsq+c2*rsq*rsq;
                    
                    return true;
                }
                else
                    return false;
            }

            DEVICE Scalar evalPressureLRCIntegral()
            {
                //throw std::runtime_error("evalPresureLRCIntegral is not implemented yet for this pair potential.");
                return 0;
            }

            DEVICE Scalar evalEnergyLRCIntegral()
            {
                //throw std::runtime_error("evalEnergyLRCIntegral is not implemented yet for this pair potential.");
                return 0;
            }

#ifndef __HIPCC__
            //! Get the name of this potential
            /*! \returns The potential name.
             */
            static std::string getName()
            {
                return std::string("polydisperse");
            }

            std::string getShapeSpec() const
            {
                throw std::runtime_error("Shape definition not supported for this pair potential.");
            }
#endif

        protected:
            Scalar rsq;     //!< Stored rsq from the constructor
            Scalar rcutsq;  //!< Stored rcutsq from the constructor
            Scalar d_i;     //!< d_i diameter of particle i
            Scalar d_j;     //!< d_j diameter of particle j
            
            //Parameters to be read
            Scalar v0;
            Scalar eps;
            Scalar rcut;

            //Additional parameters to be computed from the ones I read 
            Scalar c0;
            Scalar c1;
            Scalar c2;

            //Integer exponent of the pair potential
            int m_expnt;
            int mhalf;
            int n_expnt;
            int nhalf;
        };
    } // end namespace md
} // end namespace hoomd

#endif // __PAIR_EVALUATOR_POLYDISPERSE_H__
