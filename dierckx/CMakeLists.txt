cmake_minimum_required(VERSION 2.6)

project(dierckx)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

if(CMAKE_COMPILER_IS_GNUCC)
    set(real8_flag "-fdefault-real-8")
else()
    set(real8_flag "-r8")
endif()

set(cflags "-fPIC")

file(GLOB dierckx_SRC "*.f")

# DIERCKX
add_library(dierckx SHARED dierckx.c ${dierckx_SRC})
set_target_properties(dierckx PROPERTIES COMPILE_FLAGS ${cflags})

add_library(dierckx_static STATIC dierckx.c ${dierckx_SRC})
set_target_properties(dierckx_static PROPERTIES OUTPUT_NAME dierckx)
set_target_properties(dierckx_static PROPERTIES COMPILE_FLAGS ${cflags})

# DDIERCKX
set(dcflags "${cflags} ${real8_flag}")

add_library(ddierckx SHARED dierckx.c ${dierckx_SRC})
set_target_properties(ddierckx PROPERTIES COMPILE_FLAGS ${dcflags})

add_library(ddierckx_static STATIC dierckx.c ${dierckx_SRC})
set_target_properties(ddierckx_static PROPERTIES OUTPUT_NAME ddierckx)
set_target_properties(ddierckx_static PROPERTIES COMPILE_FLAGS ${dcflags})

add_subdirectory(python)

install(TARGETS dierckx dierckx_static ddierckx ddierckx_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES dierckx.h DESTINATION include)
