################ Python only modules
set(pyfiles __init__.py
          pair.py
   )

install(FILES ${pyfiles}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/polymd/pair
       )

MACRO(copy_file file)
    add_custom_command (
        OUTPUT ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}
        COMMENT    "Copy hoomd/polymd/pair/${file}"
    )
ENDMACRO(copy_file)
