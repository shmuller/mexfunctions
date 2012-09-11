cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT dblMaxwmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS dblMaxwmodule.c
    COMMENT "Building Python extension")

add_custom_target(dblMaxwmodule DEPENDS dblMaxwmodule_command)

add_dependencies(python dblMaxwmodule)

