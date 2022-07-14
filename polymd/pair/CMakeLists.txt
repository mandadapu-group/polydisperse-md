set(files __init__.py
          pair.py
   )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/polymd/pair
       )

copy_files_to_build("${files}" "polymd-pair" "*.py")
