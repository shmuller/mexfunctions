cmake_minimum_required(VERSION 2.6)

project(pppack)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

set(pppack_SRC pppack.f90 l2.f90)

add_library(pppack SHARED ${pppack_SRC})

add_library(pppack_static STATIC ${pppack_SRC})
set_target_properties(pppack_static PROPERTIES OUTPUT_NAME pppack)
set_target_properties(pppack_static PROPERTIES COMPILE_FLAGS "-fPIC")

add_executable(pppack_prb pppack_prb.f90)
target_link_libraries(pppack_prb pppack)

add_subdirectory(python)

install(TARGETS pppack pppack_static pppack_prb
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES pppack.h DESTINATION include)
