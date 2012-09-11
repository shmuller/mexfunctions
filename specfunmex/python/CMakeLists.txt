cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT specfunmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS specfunmodule.c
    COMMENT "Building Python extension")

add_custom_target(specfunmodule DEPENDS specfunmodule_command)

add_dependencies(python specfunmodule)

