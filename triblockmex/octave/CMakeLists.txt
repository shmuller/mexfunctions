cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT triblockmex.oct
    COMMAND mkoctfile -o triblockmex.oct
        ${CMAKE_CURRENT_SOURCE_DIR}/triblockoct.cc
        -I${_triblock_SOURCE_DIR} -L${_triblock_BINARY_DIR} -l_triblock
    DEPENDS triblockoct.cc _triblock
    COMMENT "Building Octave extension")

add_custom_target(triblockmex ALL
    DEPENDS triblockmex.oct)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/triblockmex.oct DESTINATION mex)

