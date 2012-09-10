cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT cardanomex.mex
    COMMAND mkoctfile --mex 
        ${CMAKE_CURRENT_SOURCE_DIR}/cardanomex.c
        -I${cardano_SOURCE_DIR} -L${cardano_BINARY_DIR} -lcardano
    DEPENDS cardanomex.c cardano
    COMMENT "Building Octave extension")

add_custom_target(cardanomex ALL
    DEPENDS cardanomex.mex)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/cardanomex.mex DESTINATION mex)


