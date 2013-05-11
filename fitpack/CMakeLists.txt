cmake_minimum_required(VERSION 2.6)

project(fitpack)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

add_library(fitpack SHARED fitpack.f)

add_library(fitpack_static STATIC fitpack.f)
set_target_properties(fitpack_static PROPERTIES OUTPUT_NAME fitpack)
set_target_properties(fitpack_static PROPERTIES COMPILE_FLAGS "-fPIC")

add_subdirectory(python)

install(TARGETS fitpack fitpack_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES fitpack.h DESTINATION include)
