cmake_minimum_required(VERSION 2.0)

project(dblMaxwmex)

add_executable(angle_int_test angle_int_test.c ../quadmex/gauss_legendre.c)
target_link_libraries(angle_int_test specfun common)


add_custom_command(
    OUTPUT ${dblMaxwmex_BINARY_DIR}/dblMaxwmex.mex
    COMMAND mkoctfile --mex 
        ${CMAKE_CURRENT_SOURCE_DIR}/dblMaxwmex.c 
        ${CMAKE_CURRENT_SOURCE_DIR}/../quadmex/gauss_legendre.c
        -lspecfun -latomic -lcommon
    DEPENDS dblMaxwmex.c angle_int.c
    COMMENT "Building Octave extension")

add_custom_target(dblMaxwmex ALL
    DEPENDS ${dblMaxwmex_BINARY_DIR}/dblMaxwmex.mex)

