set(COMPONENT_NAME polymd)

set(_polymd_sources module-polymd.cc)

set(_polymd_headers EvaluatorPairPolydisperse.h
    )

#set(_polymd_cu_sources 
#    AllDriverPotentialPairPluginGPU.cu
#            LJPluginDriverPotentialPairGPU.cu
#            ForceShiftedLJPluginDriverPotentialPairGPU.cu
#            PolydisperseDriverPotentialPairGPU.cu
#            PolydisperseLJDriverPotentialPairGPU.cu
#            Polydisperse18DriverPotentialPairGPU.cu
#            Polydisperse10DriverPotentialPairGPU.cu
#            PolydisperseLJ106DriverPotentialPairGPU.cu
#    )

#More Pair Potentials can be added (in theory)
set(_pair_evaluators LJ
                    #Polydisperse
                     )

foreach(_evaluator ${_pair_evaluators})
    set(_evaluator_cpp ${_evaluator})

    configure_file(export_PotentialPair.cc.inc
                   export_PotentialPair${_evaluator}.cc
                   @ONLY)
    set(_polymd_sources ${_polymd_sources} export_PotentialPair${_evaluator}.cc)

    if (ENABLE_GPU)
       configure_file(export_PotentialPairGPU.cc.inc
                       export_PotentialPair${_evaluator}GPU.cc
                       @ONLY)
       configure_file(PotentialPairGPUKernel.cu.inc
                        PotentialPair${_evaluator}GPUKernel.cu
                        @ONLY)
        set(_cuda_sources ${_cuda_sources}
            export_PotentialPair${_evaluator}GPU.cc
            PotentialPair${_evaluator}GPUKernel.cu
            )
    endif()
endforeach()

pybind11_add_module(_polymd SHARED ${_polymd_sources} ${_polymd_headers} ${_cuda_sources} NO_EXTRAS)

# alias into the HOOMD namespace so that plugins and symlinked components both work
add_library(HOOMD::_polymd ALIAS _polymd)

# specify required include directories
target_include_directories(_polymd PUBLIC ${PROJECT_SOURCE_DIR}/polymd)

if(APPLE)
set_target_properties(_polymd PROPERTIES INSTALL_RPATH "@loader_path/..;@loader_path")
else()
set_target_properties(_polymd PROPERTIES INSTALL_RPATH "\$ORIGIN/..;\$ORIGIN")
endif()

# link the library to its dependencies
if (CUSOLVER_AVAILABLE)
    # CUDA 8.0 requires that we link in gomp
    target_link_libraries(_polymd PUBLIC HOOMD::_hoomd HOOMD::_md CUDA::cusolver CUDA::cusparse gomp)
else()
    target_link_libraries(_polymd PUBLIC HOOMD::_hoomd HOOMD::_md)
endif()
if (ENABLE_GPU)
    target_link_libraries(_polymd PRIVATE HOOMD::neighbor)
endif()

fix_cudart_rpath(_polymd)

# install the library
install(TARGETS _polymd EXPORT HOOMDTargets
        LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}/polymd
        )

################ Python only modules
# copy python modules to the build directory to make it a working python package
#set(files __init__.py
#    )

#install(FILES ${files}
#        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/polymd
#       )

#copy_files_to_build("${files}" "polymd" "*.py")

# install headers in installation target
#install(FILES ${_polymd_headers}
#        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/include/hoomd/polymd
#       )

add_subdirectory(pair)
