cmake_minimum_required(VERSION 2.6)

add_custom_command(
    OUTPUT mdsclientmodule_command
    COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && rm -rf build && python setup.py install
    DEPENDS mdsclient/mdsclientmodule.c
    COMMENT "Building Python extension")

add_custom_target(mdsclientmodule DEPENDS mdsclientmodule_command)

add_dependencies(python mdsclientmodule)

