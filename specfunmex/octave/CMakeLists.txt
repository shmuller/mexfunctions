cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT specfunmex.oct
    COMMAND mkoctfile -o specfunmex.oct
        ${CMAKE_CURRENT_SOURCE_DIR}/specfunoct.cc
        -I${specfun_SOURCE_DIR} -L${specfun_BINARY_DIR} -lspecfun 
        -I${common_SOURCE_DIR} -L${common_BINARY_DIR} -lcommon
    DEPENDS specfunoct.cc specfun common
    COMMENT "Building Octave extension")

add_custom_target(specfunmex ALL
    DEPENDS specfunmex.oct)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/specfunmex.oct DESTINATION mex)

