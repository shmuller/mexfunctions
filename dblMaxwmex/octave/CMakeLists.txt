cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT dblMaxwmex.mex
    COMMAND mkoctfile --mex 
        ${CMAKE_CURRENT_SOURCE_DIR}/dblMaxwmex.c
        -I${dblMaxw_SOURCE_DIR} -L${dblMaxw_BINARY_DIR} -ldblMaxw
    DEPENDS dblMaxwmex.c dblMaxw
    COMMENT "Building Octave extension")

add_custom_target(dblMaxwmex ALL
    DEPENDS dblMaxwmex.mex)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/dblMaxwmex.mex DESTINATION mex)


