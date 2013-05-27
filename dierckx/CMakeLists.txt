cmake_minimum_required(VERSION 2.6)

project(dierckx)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

if(CMAKE_COMPILER_IS_GNUCC)
    set(real8_flag "-fdefault-real-8")
else()
    set(real8_flag "-r8")
endif()

file(GLOB dierckx_SRC "*.f")

add_library(dierckx SHARED dierckx.c ${dierckx_SRC})
set_target_properties(dierckx PROPERTIES COMPILE_FLAGS ${real8_flag})

add_library(dierckx_static STATIC dierckx.c ${dierckx_SRC})
set_target_properties(dierckx_static PROPERTIES OUTPUT_NAME dierckx)
set_target_properties(dierckx_static PROPERTIES COMPILE_FLAGS "${real8_flag} -fPIC")

add_subdirectory(python)

install(TARGETS dierckx dierckx_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES dierckx.h DESTINATION include)
