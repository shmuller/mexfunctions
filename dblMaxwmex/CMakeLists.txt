cmake_minimum_required(VERSION 2.6)

project(dblMaxw)

add_library(dblMaxw SHARED dblMaxw.c angle_int.c ../quadmex/gauss_legendre.c)
target_link_libraries(dblMaxw specfun atomic common)

add_executable(angle_int_test angle_int_test.c angle_int.c ../quadmex/gauss_legendre.c)
target_link_libraries(angle_int_test specfun common)


add_custom_command(
    OUTPUT ${dblMaxw_BINARY_DIR}/dblMaxwmex.mex
    COMMAND mkoctfile --mex 
        ${CMAKE_CURRENT_SOURCE_DIR}/dblMaxwmex.c
        ${CMAKE_CURRENT_SOURCE_DIR}/dblMaxw.c
        ${CMAKE_CURRENT_SOURCE_DIR}/angle_int.c
        ${CMAKE_CURRENT_SOURCE_DIR}/../quadmex/gauss_legendre.c
        -lspecfun -latomic -lcommon
    DEPENDS dblMaxwmex.c dblMaxw.c angle_int.c
    COMMENT "Building Octave extension")

#add_custom_command(
#    OUTPUT ${dblMaxw_BINARY_DIR}/dblMaxwmex.mex
#    COMMAND mkoctfile --mex 
#        ${CMAKE_CURRENT_SOURCE_DIR}/dblMaxwmex.c 
#        -L${dblMaxw_BINARY_DIR} -ldblMaxw
#    DEPENDS dblMaxwmex.c dblMaxw
#    COMMENT "Building Octave extension")

add_custom_target(dblMaxwmex ALL
    DEPENDS ${dblMaxw_BINARY_DIR}/dblMaxwmex.mex)


set(CMAKE_BUILD_TYPE Release)

