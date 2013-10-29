cmake_minimum_required(VERSION 2.8.8)

project(slatec)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

if(CMAKE_COMPILER_IS_GNUCC)
    set(real8_flag "-fdefault-real-8")
else()
    set(real8_flag "-r8")
endif()

#set(cflags "${real8_flag} -fPIC")
set(cflags "-fPIC")
set(lflags "-undefined dynamic_lookup")

file(GLOB slatec_SRC "src/*.f")
set(mach_SRC d1mach.f90 i1mach.f90 r1mach.f90)

add_library(slatec_obj OBJECT ${slatec_SRC} ${mach_SRC})
set_target_properties(slatec_obj PROPERTIES COMPILE_FLAGS ${cflags})

add_library(slatec_static STATIC $<TARGET_OBJECTS:slatec_obj>)
set_target_properties(slatec_static PROPERTIES OUTPUT_NAME slatec)

add_library(slatec SHARED $<TARGET_OBJECTS:slatec_obj>)
set_target_properties(slatec PROPERTIES LINK_FLAGS ${lflags})

file(GLOB slatec_chk_SRC "chk/*.f")
add_library(slatec_chk SHARED ${slatec_chk_SRC})
set_target_properties(slatec_chk PROPERTIES COMPILE_FLAGS ${cflags})
target_link_libraries(slatec_chk slatec)

file(GLOB test_src "chk/tests/*.f")
foreach(testsourcefile ${test_src})
    get_filename_component(testname ${testsourcefile} NAME_WE)
    message(STATUS "Adding executable ${testname} from source ${testsourcefile}")
    add_executable(${testname} ${testsourcefile})
    set_target_properties(${testname} PROPERTIES COMPILE_FLAGS ${cflags})
    target_link_libraries(${testname} slatec_chk)
endforeach(testsourcefile)

install(TARGETS slatec slatec_static slatec_chk
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES slatec.h DESTINATION include)
