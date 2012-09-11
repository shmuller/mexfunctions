cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT mediansmoothmex.mex
    COMMAND mkoctfile --mex 
        ${CMAKE_CURRENT_SOURCE_DIR}/mediansmoothmex.c
    DEPENDS mediansmoothmex.c ../mediansmooth.c
    COMMENT "Building Octave extension")

add_custom_target(mediansmoothmex ALL
    DEPENDS mediansmoothmex.mex)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mediansmoothmex.mex DESTINATION mex)

install(FILES mediansmooth.m DESTINATION mex)

