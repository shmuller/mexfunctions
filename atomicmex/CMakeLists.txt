cmake_minimum_required(VERSION 2.0)

project(atomic)

add_library(atomic SHARED atomic.c)
target_link_libraries(atomic common)

add_library(atomic_static STATIC atomic.c)
set_target_properties(atomic_static PROPERTIES OUTPUT_NAME atomic)
set_target_properties(atomic_static PROPERTIES COMPILE_FLAGS -fPIC)

install(TARGETS atomic atomic_static 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

