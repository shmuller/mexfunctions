cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT triblockmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS triblockmodule.c
    COMMENT "Building Python extension")

add_custom_target(triblockmodule DEPENDS triblockmodule_command)

add_dependencies(python triblockmodule)

