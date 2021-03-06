cmake_minimum_required(VERSION 2.8.8)

project(slatec)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_C_FLAGS_RELEASE "-O4 -finline-functions")

enable_language(Fortran)

if(CMAKE_COMPILER_IS_GNUCC)
    set(real8_flag "-fdefault-real-8")
else()
    set(real8_flag "-r8")
endif()

set(cflags "-fPIC")
set(lflags "-undefined dynamic_lookup")

file(GLOB slatec_SRC "src/*.f")
file(GLOB extra_SRC "ext/*.f90")
set(mach_SRC d1mach.f90 i1mach.f90 r1mach.f90)

# SLATEC
add_library(slatec_obj OBJECT ${slatec_SRC} ${mach_SRC})
set_target_properties(slatec_obj PROPERTIES COMPILE_FLAGS ${cflags})

add_library(slatec_ext OBJECT ${extra_SRC})
set_target_properties(slatec_ext PROPERTIES COMPILE_FLAGS "${cflags} -fimplicit-none")

add_library(slatec_static STATIC $<TARGET_OBJECTS:slatec_obj> $<TARGET_OBJECTS:slatec_ext>)
set_target_properties(slatec_static PROPERTIES OUTPUT_NAME slatec)

add_library(slatec SHARED $<TARGET_OBJECTS:slatec_obj> $<TARGET_OBJECTS:slatec_ext>)
#set_target_properties(slatec PROPERTIES LINK_FLAGS ${lflags})

# DSLATEC
set(dcflags "${cflags} ${real8_flag}")

add_library(dslatec_obj OBJECT ${slatec_SRC} ${mach_SRC})
set_target_properties(dslatec_obj PROPERTIES COMPILE_FLAGS ${dcflags})

add_library(dslatec_ext OBJECT ${extra_SRC})
set_target_properties(dslatec_ext PROPERTIES COMPILE_FLAGS "${dcflags} -fimplicit-none")

add_library(dslatec_static STATIC $<TARGET_OBJECTS:dslatec_obj> $<TARGET_OBJECTS:dslatec_ext>)
set_target_properties(dslatec_static PROPERTIES OUTPUT_NAME dslatec)

add_library(dslatec SHARED $<TARGET_OBJECTS:dslatec_obj> $<TARGET_OBJECTS:dslatec_ext>)
#set_target_properties(dslatec PROPERTIES LINK_FLAGS ${lflags})

# SLATEC TESTS
file(GLOB slatec_chk_SRC "chk/*.f")
add_library(slatec_chk SHARED ${slatec_chk_SRC})
set_target_properties(slatec_chk PROPERTIES COMPILE_FLAGS ${cflags})
#set_target_properties(slatec_chk PROPERTIES LINK_FLAGS ${lflags})
target_link_libraries(slatec_chk slatec)

file(GLOB test_src "chk/tests/*.f")
foreach(testsourcefile ${test_src})
    get_filename_component(testname ${testsourcefile} NAME_WE)
    message(STATUS "Adding executable ${testname} from source ${testsourcefile}")
    add_executable(${testname} ${testsourcefile})
    set_target_properties(${testname} PROPERTIES COMPILE_FLAGS ${cflags})
    target_link_libraries(${testname} slatec_chk)
endforeach(testsourcefile)

install(TARGETS slatec slatec_static dslatec dslatec_static slatec_chk
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES slatec.h DESTINATION include)
