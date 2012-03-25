cmake_minimum_required(VERSION 2.0)

project(angle_int_test)

add_executable(angle_int_test angle_int_test.c ../quadmex/gauss_legendre.c)
target_link_libraries(angle_int_test specfun common)

