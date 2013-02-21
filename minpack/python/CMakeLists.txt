cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT minpackmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS minpackmodule.c
    COMMENT "Building Python extension")

add_custom_target(minpackmodule DEPENDS minpackmodule_command)

add_dependencies(python minpackmodule)

