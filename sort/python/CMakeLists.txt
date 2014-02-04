cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT sortmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS sortmodule.c
    COMMENT "Building Python extension")

add_custom_target(sortmodule DEPENDS sortmodule_command)

add_dependencies(python sortmodule)

