cmake_minimum_required(VERSION 2.6)

project(atomic)

set(CMAKE_BUILD_TYPE Release)

add_library(atomic SHARED atomic.c)
target_link_libraries(atomic common)

add_library(atomic_static STATIC atomic.c)
set_target_properties(atomic_static PROPERTIES OUTPUT_NAME atomic)
set_target_properties(atomic_static PROPERTIES COMPILE_FLAGS -fPIC)

set(CMAKE_BUILD_TYPE Release)

install(TARGETS atomic atomic_static 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

