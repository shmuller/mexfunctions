cmake_minimum_required(VERSION 2.6)

project(mexfunctions)

add_custom_target(python)

add_subdirectory(common)
add_subdirectory(specfunmex)
add_subdirectory(atomicmex)
add_subdirectory(dblMaxwmex)
add_subdirectory(cardanomex)
add_subdirectory(mdsclient)
add_subdirectory(mediansmoothmex)

set(CMAKE_BUILD_TYPE Release)

