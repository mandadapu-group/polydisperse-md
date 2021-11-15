set(COMPONENT_NAME polymd)

set(_${COMPONENT_NAME}_sources 
                    module-md-plugin.cc
                    )

set(_${COMPONENT_NAME}_cu_sources 
            LJPluginDriverPotentialPairGPU.cu
            ForceShiftedLJPluginDriverPotentialPairGPU.cu
            PolydisperseDriverPotentialPairGPU.cu
            PolydisperseLJDriverPotentialPairGPU.cu
            Polydisperse18DriverPotentialPairGPU.cu
            Polydisperse10DriverPotentialPairGPU.cu
            PolydisperseLJ106DriverPotentialPairGPU.cu
    )

# Need to define NO_IMPORT_ARRAY in every file but module.cc
set_source_files_properties(${_${COMPONENT_NAME}_sources} ${_${COMPONENT_NAME}_cu_sources} PROPERTIES COMPILE_DEFINITIONS NO_IMPORT_ARRAY)
                
if (ENABLE_CUDA)
CUDA_COMPILE(_CUDA_GENERATED_FILES ${_${COMPONENT_NAME}_cu_sources} OPTIONS ${CUDA_ADDITIONAL_OPTIONS} SHARED)
endif (ENABLE_CUDA)

pybind11_add_module (_${COMPONENT_NAME} SHARED ${_${COMPONENT_NAME}_sources} ${_CUDA_GENERATED_FILES} NO_EXTRAS)
#pybind11_add_module (_${COMPONENT_NAME} SHARED ${_${COMPONENT_NAME}_sources} NO_EXTRAS)

# link the library to its dependencies
target_link_libraries(_${COMPONENT_NAME} ${HOOMD_LIBRARIES} ${HOOMD_MD_LIB})

# if we are compiling with MPI support built in, set appropriate
# compiler/linker flags
if (ENABLE_MPI)
   if(MPI_COMPILE_FLAGS)
       set_target_properties(_${COMPONENT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
   endif(MPI_COMPILE_FLAGS)
   if(MPI_LINK_FLAGS)
       set_target_properties(_${COMPONENT_NAME} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
   endif(MPI_LINK_FLAGS)
endif(ENABLE_MPI)

if(APPLE)
    set_target_properties(_${COMPONENT_NAME} PROPERTIES INSTALL_RPATH "@loader_path/..;@loader_path")
    set_target_properties(_${COMPONENT_NAME} PROPERTIES SUFFIX ".so")
else()
    set_target_properties(_${COMPONENT_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/..;$ORIGIN")
endif()

fix_cudart_rpath(_${COMPONENT_NAME})

# install the library
install(TARGETS _${COMPONENT_NAME} 
        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
        )

################ Python only modules
# copy python modules to the build directory to make it a working python package
MACRO(copy_file file)
    add_custom_command (
        OUTPUT ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}
        COMMENT    "Copy hoomd/${COMPONENT_NAME}/${file}"
    )
ENDMACRO(copy_file)

set(files   __init__.py
            pair.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/${COMPONENT_NAME}
       )

foreach(file ${files})
    copy_file(${file})
endforeach()

add_custom_target(copy_${COMPONENT_NAME} ALL DEPENDS ${files})
