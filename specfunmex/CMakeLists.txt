cmake_minimum_required(VERSION 2.6)

project(specfun)
enable_language(Fortran)

add_library(specfun SHARED specfun.f90)
set_target_properties(specfun PROPERTIES COMPILE_FLAGS "-fno-underscoring")
target_link_libraries(specfun common)

add_library(specfun_static STATIC specfun.f90)
set_target_properties(specfun_static PROPERTIES OUTPUT_NAME specfun)
set_target_properties(specfun_static PROPERTIES COMPILE_FLAGS "-fPIC -fno-underscoring")

add_subdirectory(octave)
add_subdirectory(python)

install(TARGETS specfun specfun_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

