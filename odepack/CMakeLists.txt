cmake_minimum_required(VERSION 2.6)

project(odepack)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

set(odepack_SRC opkdmain.f opkda1.f opkda2.f)

add_library(odepack SHARED ${odepack_SRC})

add_library(odepack_static STATIC ${odepack_SRC})
set_target_properties(odepack_static PROPERTIES OUTPUT_NAME odepack)
set_target_properties(odepack_static PROPERTIES COMPILE_FLAGS "-fPIC")

add_subdirectory(python)

install(TARGETS odepack odepack_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES odepack.h DESTINATION include)
