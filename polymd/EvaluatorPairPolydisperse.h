// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander

#ifndef __PAIR_EVALUATOR_POLYDISPERSE_H__
#define __PAIR_EVALUATOR_POLYDISPERSE_H__

#ifndef NVCC
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairPolydisperse.h
    \brief Defines the pair evaluator class for LJ potentials
    \details As the prototypical example of a MD pair potential, this also serves as the primary documentation and
    base reference for the implementation of pair evaluators.
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

//! Class for evaluating the LJ pair potential
/*! <b>General Overview</b>

    EvaluatorPairPolydisperse is a low level computation class that computes the LJ pair potential V(r). As the standard
    MD potential, it also serves as a well documented example of how to write additional pair potentials. "Standard"
    pair potentials in hoomd are all handled via the template class PotentialPair. PotentialPair takes a potential
    evaluator as a template argument. In this way, all the complicated data management and other details of computing
    the pair force and potential on every single particle is only written once in the template class and the difference
    V(r) potentials that can be calculated are simply handled with various evaluator classes. Template instantiation is
    equivalent to inlining code, so there is no performance loss.

    In hoomd, a "standard" pair potential is defined as V(rsq, rcutsq, params, di, dj, qi, qj), where rsq is the squared
    distance between the two particles, rcutsq is the cutoff radius at which the potential goes to 0, params is any
    number of per type-pair parameters, di, dj are the diameters of particles i and j, and qi, qj are the charges of
    particles i and j respectively.

    Diameter and charge are not always needed by a given pair evaluator, so it must provide the functions
    needsDiameter() and needsCharge() which return boolean values signifying if they need those quantities or not. A
    false return value notifies PotentialPair that it need not even load those values from memory, boosting performance.

    If needsDiameter() returns true, a setDiameter(Scalar di, Scalar dj) method will be called to set the two diameters.
    Similarly, if needsCharge() returns true, a setCharge(Scalar qi, Scalar qj) method will be called to set the two
    charges.

    All other arguments are common among all pair potentials and passed into the constructor. Coefficients are handled
    in a special way: the pair evaluator class (and PotentialPair) manage only a single parameter variable for each
    type pair. Pair potentials that need more than 1 parameter can specify that their param_type be a compound
    structure and reference that. For coalesced read performance on G200 GPUs, it is highly recommended that param_type
    is one of the following types: Scalar, Scalar2, Scalar4.

    The program flow will proceed like this: When a potential between a pair of particles is to be evaluated, a
    PairEvaluator is instantiated, passing the common parameters to the constructor and calling setDiameter() and/or
    setCharge() if need be. Then, the evalForceAndEnergy() method is called to evaluate the force and energy (more
    on that later). Thus, the evaluator must save all of the values it needs to compute the force and energy in member
    variables.

    evalForceAndEnergy() makes the necessary computations and sets the out parameters with the computed values.
    Specifically after the method completes, \a force_divr must be set to the value
    \f$ -\frac{1}{r}\frac{\partial V}{\partial r}\f$ and \a pair_eng must be set to the value \f$ V(r) \f$ if \a energy_shift is false or
    \f$ V(r) - V(r_{\mathrm{cut}}) \f$ if \a energy_shift is true.

    A pair potential evaluator class is also used on the GPU. So all of its members must be declared with the
    DEVICE keyword before them to mark them __device__ when compiling in nvcc and blank otherwise. If any other code
    needs to diverge between the host and device (i.e., to use a special math function like __powf on the device), it
    can similarly be put inside an ifdef NVCC block.

    <b>LJ specifics</b>

    EvaluatorPairPolydisperse evaluates the function:
    \f[ V_{\mathrm{LJ}}(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
                                            \alpha \left( \frac{\sigma}{r} \right)^{6} \right] \f]
    broken up as follows for efficiency
    \f[ V_{\mathrm{LJ}}(r) = r^{-6} \cdot \left( 4 \varepsilon \sigma^{12} \cdot r^{-6} -
                                            4 \alpha \varepsilon \sigma^{6} \right) \f]
    . Similarly,
    \f[ -\frac{1}{r} \frac{\partial V_{\mathrm{LJ}}}{\partial r} = r^{-2} \cdot r^{-6} \cdot
            \left( 12 \cdot 4 \varepsilon \sigma^{12} \cdot r^{-6} - 6 \cdot 4 \alpha \varepsilon \sigma^{6} \right) \f]

    The LJ potential does not need diameter or charge. Two parameters are specified and stored in a Scalar2. \a lj1 is
    placed in \a params.x and \a lj2 is in \a params.y.

    These are related to the standard lj parameters sigma and epsilon by:
    - \a lj1 = 4.0 * epsilon * pow(sigma,12.0)
    - \a lj2 = alpha * 4.0 * epsilon * pow(sigma,6.0);

*/
class EvaluatorPairPolydisperse
    {
    public:
        //! Define the parameter type used by this pair potential evaluator
        typedef Scalar3 param_type;

        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        DEVICE EvaluatorPairPolydisperse(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
            : rsq(_rsq), rcutsq(_rcutsq), v0(_params.x), eps(_params.y), scaledr_cut(_params.z), c0(Scalar(-1.92415)), c1(Scalar(2.11106)), c2(Scalar(-0.591097))
            {
                c0 =  Scalar(-28.0)*v0/pow(scaledr_cut,12);
                c1 =  Scalar(48.0)*v0/pow(scaledr_cut,14);
                c2 =  Scalar(-21.0)*v0/pow(scaledr_cut,16);
            }

        //! LJ doesn't use diameter
        DEVICE static bool needsDiameter() { return true; }
        //! Accept the optional diameter values
        /*! \param di Diameter of particle i
            \param dj Diameter of particle j
        */
        DEVICE void setDiameter(Scalar di, Scalar dj) 
        {
           d_i = di;
           d_j = dj; 
        }

        //! LJ doesn't use charge
        DEVICE static bool needsCharge() { return false; }
        //! Accept the optional diameter values
        /*! \param qi Charge of particle i
            \param qj Charge of particle j
        */
        DEVICE void setCharge(Scalar qi, Scalar qj) { }

        //! Evaluate the force and energy
        /*! \param force_divr Output parameter to write the computed force divided by r.
            \param pair_eng Output parameter to write the computed pair energy
            \param energy_shift If true, the potential must be shifted so that V(r) is continuous at the cutoff
            \note There is no need to check if rsq < rcutsq in this method. Cutoff tests are performed
                  in PotentialPair.

            \return True if they are evaluated or false if they are not because we are beyond the cutoff
        */
        DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
            {
                Scalar sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
                Scalar actualcutsq = scaledr_cut*scaledr_cut*sigma*sigma;
                // compute the force divided by r in force_divr
                if (rsq < actualcutsq && v0 != 0)
                    {
                    Scalar r2inv = sigma*sigma*Scalar(1.0)/rsq;
                    Scalar _rsq = Scalar(1.0)*rsq/(sigma*sigma);
                    Scalar r6inv = r2inv * r2inv * r2inv;
                    force_divr = (Scalar(12.0)*v0*r2inv*r6inv*r6inv-Scalar(2.0)*c1 -Scalar(4.0)*c2*_rsq)/(sigma*sigma);
                    
                    //No energy shift is needed
                    pair_eng = v0*r6inv*r6inv+c0+c1*_rsq+c2*_rsq*_rsq;
                    
                    return true;
                    }
                else
                    return false;
            }

        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Must be short and all lowercase, as this is the name energies will be logged as
            via analyze.log.
        */
        static std::string getName()
            {
            return std::string("polydisperse-12");
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
        Scalar scaledr_cut;

        //Additional parameters to be computed from the ones I read 
        Scalar c0;
        Scalar c1;
        Scalar c2;
    };

#endif // __PAIR_EVALUATOR_POLYDISPERSE_H__
