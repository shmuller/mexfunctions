cmake_minimum_required(VERSION 2.6)

project(sort)

set(CMAKE_BUILD_TYPE Release)

add_library(sort SHARED sort.c sort.cpp)

add_subdirectory(python)

set(CMAKE_BUILD_TYPE Release)

install(TARGETS sort
    RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(FILES sort.h DESTINATION include)

