cmake_minimum_required(VERSION 2.6)

project(fitfun)

set(CMAKE_BUILD_TYPE Release)

add_library(fitfun SHARED fitfun.c)
target_link_libraries(fitfun minpack)

add_subdirectory(python)

install(TARGETS fitfun
    RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(FILES fitfun.h DESTINATION include)

