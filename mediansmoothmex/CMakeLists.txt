cmake_minimum_required(VERSION 2.6)

project(mediansmooth)

add_library(mediansmooth SHARED mediansmooth.c mediansmooth_pqueue.c pqueue2.c pqueue.c)
add_dependencies(mediansmooth mediansmooth_template.c)

add_subdirectory(octave)
add_subdirectory(python)

set(CMAKE_BUILD_TYPE Release)

install(TARGETS mediansmooth
    RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(FILES mediansmooth.h DESTINATION include)

