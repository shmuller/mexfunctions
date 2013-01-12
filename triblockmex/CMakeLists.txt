cmake_minimum_required(VERSION 2.6)

project(_triblock)

set(CMAKE_BUILD_TYPE Release)

add_library(_triblock SHARED _triblock.c)
target_link_libraries(_triblock lapack blas)

add_library(_triblock_static STATIC _triblock.c)
set_target_properties(_triblock_static PROPERTIES OUTPUT_NAME _triblock)
set_target_properties(_triblock_static PROPERTIES COMPILE_FLAGS -fPIC)

if(mkoctfile_EXEC)
add_subdirectory(octave)
endif(mkoctfile_EXEC)

add_subdirectory(python)

install(TARGETS _triblock _triblock_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES _triblock.h DESTINATION include)

