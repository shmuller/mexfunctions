cmake_minimum_required(VERSION 2.6)

project(minpack)
enable_language(Fortran)

set(CMAKE_BUILD_TYPE Release)

add_library(minpack SHARED minpack.c dpmpar.f enorm.f fdjac2.f lmder.f 
    lmdif.f lmdif1.f lmpar.f qrfac.f qrsolv.f)

add_library(minpack_static STATIC minpack.c dpmpar.f enorm.f fdjac2.f lmder.f 
    lmdif.f lmdif1.f lmpar.f qrfac.f qrsolv.f)
set_target_properties(minpack_static PROPERTIES OUTPUT_NAME minpack)
set_target_properties(minpack_static PROPERTIES COMPILE_FLAGS "-fPIC")

add_executable(minpack_test minpack_test.c)
target_link_libraries(minpack_test minpack)

add_subdirectory(python)

install(TARGETS minpack minpack_static minpack_test
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES minpack.h DESTINATION include)
