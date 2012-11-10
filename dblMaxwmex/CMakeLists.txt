cmake_minimum_required(VERSION 2.6)

project(dblMaxw)

add_library(dblMaxw SHARED dblMaxw.c angle_int.c ../quadmex/gauss_legendre.c)
target_link_libraries(dblMaxw specfun atomic common)

add_executable(angle_int_test angle_int_test.c)
target_link_libraries(angle_int_test dblMaxw)

if(mkoctfile_EXEC)
add_subdirectory(octave)
endif(mkoctfile_EXEC)

add_subdirectory(python)

install(TARGETS dblMaxw angle_int_test
    RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(FILES dblMaxw.h angle_int.h DESTINATION include)

