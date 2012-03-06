cmake_minimum_required(VERSION 2.0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")

project(mdsclient)

set(mdsplus_USE_STATIC "ON" CACHE BOOL "Link against static mdsplus libraries.")

find_package(mdsplus REQUIRED)

include_directories(${mdsplus_INCLUDE_DIRS})

add_library(mdsclient SHARED mdsclient.c tcp.c)
target_link_libraries(mdsclient ${mdsplus_LIBRARIES})

add_library(mdsclient_static STATIC mdsclient.c tcp.c)
set_target_properties(mdsclient_static PROPERTIES OUTPUT_NAME mdsclient)
set_target_properties(mdsclient_static PROPERTIES COMPILE_FLAGS -fPIC)
target_link_libraries(mdsclient_static ${mdsplus_LIBRARIES})

add_executable(mdsclient_test mdsclient_test.c)
target_link_libraries(mdsclient_test mdsclient)

install(TARGETS mdsclient_test mdsclient mdsclient_static 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

