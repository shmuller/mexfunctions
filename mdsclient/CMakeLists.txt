cmake_minimum_required(VERSION 2.0)

project(mdsclient)

set(CMAKE_BUILD_TYPE Release)

find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

add_library(mdsclient SHARED mdsclient.c tcp.c mdsiputil.c)
target_link_libraries(mdsclient ${ZLIB_LIBRARIES})

add_library(mdsclient_static STATIC mdsclient.c tcp.c mdsiputil.c)
set_target_properties(mdsclient_static PROPERTIES OUTPUT_NAME mdsclient)
set_target_properties(mdsclient_static PROPERTIES COMPILE_FLAGS -fPIC)
target_link_libraries(mdsclient_static ${ZLIB_LIBRARIES})

add_executable(mdsclient_test mdsclient_test.c)
target_link_libraries(mdsclient_test mdsclient)
    
if(mkoctfile_EXEC)
add_subdirectory(octave)
endif(mkoctfile_EXEC)

add_subdirectory(python)

install(TARGETS mdsclient_test mdsclient mdsclient_static 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(FILES mdsclient.h DESTINATION include)
