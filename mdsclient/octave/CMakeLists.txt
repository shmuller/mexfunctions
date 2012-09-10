cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT mdsclientmex.oct
    COMMAND mkoctfile -o mdsclientmex.oct
        ${CMAKE_CURRENT_SOURCE_DIR}/mdsclientoct.cc
        -I${mdsclient_SOURCE_DIR} -L${mdsclient_BINARY_DIR} -lmdsclient
    DEPENDS mdsclientoct.cc mdsclient
    COMMENT "Building Octave extension")

add_custom_target(mdsclientmex ALL
    DEPENDS mdsclientmex.oct)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mdsclientmex.oct DESTINATION mex)


