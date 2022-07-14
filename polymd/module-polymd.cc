// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace hoomd
{
    namespace md
    {
        namespace detail
        {
        //! Create the python module
        /*! each class setup their own python exports in a function export_ClassName
            create the hoomd python module and define the exports here.
        */
            void export_PotentialPairPolydisperse(pybind11::module& m);
            void export_PotentialPairLJ(pybind11::module& m);

#ifdef ENABLE_GPU
            void export_PotentialPairPolydisperse(pybind11::module& m);
            void export_PotentialPairLJGPU(pybind11::module& m);
#endif
        } // namespace detail
    } // namespace polymd
} // namespace hoomd



using namespace hoomd;
using namespace hoomd::md;
using namespace hoomd::md::detail;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
    create the md python module and define the exports here.
*/
PYBIND11_MODULE(_polymd, m)
{
    export_PotentialPairPolydisperse(m);
    export_PotentialPairLJ(m);

#ifdef ENABLE_GPU
    export_PotentialPairPolydisperseGPU(m);
    export_PotentialPairLJGPU(m);
#endif

}
