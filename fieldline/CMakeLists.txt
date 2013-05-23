cmake_minimum_required(VERSION 2.6)

project(fieldline)

set(CMAKE_BUILD_TYPE Release)

set(fieldline_SRC fieldline.c)

add_library(fieldline SHARED ${fieldline_SRC})
target_link_libraries(fieldline odepack dierckx)

add_library(fieldline_static STATIC ${fieldline_SRC})
set_target_properties(fieldline_static PROPERTIES OUTPUT_NAME fieldline)
set_target_properties(fieldline_static PROPERTIES COMPILE_FLAGS "-fPIC")
target_link_libraries(fieldline_static odepack dierckx)

add_subdirectory(python)

install(TARGETS fieldline fieldline_static
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES fieldline.h DESTINATION include)
