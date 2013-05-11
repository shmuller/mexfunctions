cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT fitpackmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS
    COMMENT "Building Python extension")

add_custom_target(fitpackmodule DEPENDS fitpackmodule_command)

add_dependencies(python fitpackmodule)

