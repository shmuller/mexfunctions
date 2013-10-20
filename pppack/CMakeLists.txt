cmake_minimum_required(VERSION 2.6)

project(pppack)

set(CMAKE_BUILD_TYPE Release)

enable_language(Fortran)

if(CMAKE_COMPILER_IS_GNUCC)
    set(real8_flag "-fdefault-real-8")
else()
    set(real8_flag "-r8")
endif()

set(pppack_SRC
    netlib/banfac.f
    netlib/banslv.f
    netlib/bchfac.f
    netlib/bchslv.f
    netlib/bsplpp.f
    netlib/bsplvb.f
    netlib/bsplvd.f
    netlib/bspp2d.f
    netlib/bvalue.f
    netlib/chol1d.f
    netlib/colloc.f
    netlib/colpnt.f
    netlib/cubspl.f
    netlib/cwidth.f
    netlib/difequ.f
    netlib/eqblok.f
    #netlib/iiex.f
    netlib/interv.f
    #netlib/ivex.f
    #netlib/ixex.f
    netlib/knots.f
    netlib/l2appr.f
    netlib/l2err.f
    netlib/l2knts.f
    #netlib/l2main.f
    netlib/newnot.f
    #netlib/newnotfake.f
    netlib/ppvalu.f
    netlib/putit.f
    netlib/round.f
    netlib/setdatex2.f
    netlib/setdatex3.f
    netlib/setdatex4.f
    netlib/setupq.f
    netlib/slvblktext.f
    netlib/smooth.f
    netlib/spli2d.f
    netlib/splint.f
    netlib/splopt.f
    netlib/tautsp.f
    netlib/titand.f
    #netlib/xex1.f
    #netlib/xex2.f
    #netlib/xex3.f
    #netlib/xex4.f
    #netlib/xiiex2.f
    #netlib/xiiex3.f
    #netlib/xiiex4.f
    #netlib/xiiiex1.f
    #netlib/xiiiex2.f
    #netlib/xiiiex2m.f
    #netlib/xiiiex3.f
    #netlib/xivex1.f
    #netlib/xvex.f
    #netlib/xviex2.f
    #netlib/xviex3.f
    #netlib/xviiex2.f
    #netlib/xviiex3.f
    bsplppd.f
    ppual.f
    ppualder.f
)

add_library(pppack SHARED ${pppack_SRC})
set_target_properties(pppack PROPERTIES COMPILE_FLAGS ${real8_flag})

add_library(pppack_static STATIC ${pppack_SRC})
set_target_properties(pppack_static PROPERTIES OUTPUT_NAME pppack)
set_target_properties(pppack_static PROPERTIES COMPILE_FLAGS "${real8_flag} -fPIC")

add_executable(l2main netlib/l2main.f)
set_target_properties(l2main PROPERTIES COMPILE_FLAGS ${real8_flag})
target_link_libraries(l2main pppack)

add_subdirectory(python)

install(TARGETS pppack pppack_static l2main
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES pppack.h DESTINATION include)

