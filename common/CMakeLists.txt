cmake_minimum_required(VERSION 2.6)

project(common)

add_library(common SHARED common.c)

add_library(common_static STATIC common.c)
set_target_properties(common_static PROPERTIES OUTPUT_NAME common)
set_target_properties(common_static PROPERTIES COMPILE_FLAGS -fPIC)

set(CMAKE_BUILD_TYPE Release)

install(TARGETS common common_static 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

