cmake_minimum_required(VERSION 2.6)

project(cardano)

add_library(cardano SHARED cardano.c)

add_subdirectory(octave)

set(CMAKE_BUILD_TYPE Release)

install(TARGETS cardano
    RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(FILES cardano.h DESTINATION include)

