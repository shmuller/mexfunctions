cmake_minimum_required(VERSION 2.6)

project(mexfunctions)

find_program(mkoctfile_EXEC NAMES mkoctfile PATHS ENV PATH)

message(STATUS "mkoctfile_EXEC=${mkoctfile_EXEC}")

# user must explicitly type "make octave" and "make python"
add_custom_target(octave)
add_custom_target(python)

add_subdirectory(common)
add_subdirectory(specfunmex)
add_subdirectory(atomicmex)
add_subdirectory(dblMaxwmex)
add_subdirectory(cardanomex)
add_subdirectory(mdsclient)
add_subdirectory(mediansmoothmex)
add_subdirectory(triblockmex)
add_subdirectory(minpack)
add_subdirectory(fitfun)

set(CMAKE_BUILD_TYPE Release)

