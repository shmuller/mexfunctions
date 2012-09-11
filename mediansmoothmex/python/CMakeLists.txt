cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT mediansmoothmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS mediansmoothmodule.c
    COMMENT "Building Python extension")

add_custom_target(mediansmoothmodule DEPENDS mediansmoothmodule_command)

add_dependencies(python mediansmoothmodule)

