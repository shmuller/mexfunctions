cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT fitfunmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS fitfunmodule.c
    COMMENT "Building Python extension")

add_custom_target(fitfunmodule DEPENDS fitfunmodule_command)

add_dependencies(python fitfunmodule)

