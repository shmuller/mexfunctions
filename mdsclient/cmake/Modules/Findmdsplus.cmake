# - Try to find mdsplus
# Once done, this will define
#
#  mdsplus_FOUND - system has mdsplus
#  mdsplus_INCLUDE_DIRS - the mdsplus include directories
#  mdsplus_LIBRARIES - link these to use mdsplus

set(paths_inc $ENV{MDSPLUS_DIR}/include /usr/local/include /usr/local/mdsplus/include)
set(paths_lib $ENV{MDSPLUS_DIR}/lib /usr/local/lib /usr/local/mdsplus/lib)

find_path(mdsplus_INCLUDE_DIR ipdesc.h PATHS ${paths_inc})
set(mdsplus_INCLUDE_DIRS ${mdsplus_INCLUDE_DIR})

if(mdsplus_USE_STATIC)
  find_library(mdsplus_libMdsIpShr_a NAMES libMdsIpShr.a PATHS ${paths_lib})
  find_library(mdsplus_libMdsIpUtil_a NAMES libMdsIpUtil.a PATHS ${paths_lib})
  set(mdsplus_LIBRARIES ${mdsplus_libMdsIpShr_a} ${mdsplus_libMdsIpUtil_a})
else(mdsplus_USE_STATIC)
  find_library(mdsplus_LIBRARY NAMES MdsIpShr PATHS ${paths_lib})
  set(mdsplus_LIBRARIES ${mdsplus_LIBRARY})
endif(mdsplus_USE_STATIC)

