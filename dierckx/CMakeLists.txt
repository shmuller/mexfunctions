cmake_minimum_required(VERSION 2.6)

project(dierckx)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

file(GLOB dierckx_SRC "*.f")

add_library(dierckx SHARED dierckx.c ${dierckx_SRC})
set_target_properties(dierckx PROPERTIES COMPILE_FLAGS "-fdefault-real-8")

add_library(dierckx_static STATIC dierckx.c ${dierckx_SRC})
set_target_properties(dierckx_static PROPERTIES OUTPUT_NAME dierckx)
set_target_properties(dierckx_static PROPERTIES COMPILE_FLAGS "-fPIC -fdefault-real-8")

add_subdirectory(python)

install(TARGETS dierckx dierckx_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES dierckx.h DESTINATION include)
