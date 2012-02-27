# - Try to find mdsplus
# Once done, this will define
#
#  mdsplus_FOUND - system has mdsplus
#  mdsplus_INCLUDE_DIRS - the mdsplus include directories
#  mdsplus_LIBRARIES - link these to use mdsplus

set(hints_inc /usr/local/include /usr/local/mdsplus/include)
set(hints_lib /usr/local/lib /usr/local/mdsplus/lib)

find_path(mdsplus_INCLUDE_DIR ipdesc.h HINTS ${hints_inc})
set(mdsplus_INCLUDE_DIRS ${mdsplus_INCLUDE_DIR})

if(mdsplus_USE_STATIC)
  find_library(mdsplus_libMdsIpShr_a NAMES libMdsIpShr.a HINTS ${hints_lib})
  find_library(mdsplus_libMdsIpUtil_a NAMES libMdsIpUtil.a HINTS ${hints_lib})
  set(mdsplus_LIBRARIES ${mdsplus_libMdsIpShr_a} ${mdsplus_libMdsIpUtil_a})
else(mdsplus_USE_STATIC)
  find_library(mdsplus_LIBRARY NAMES MdsIpShr HINTS ${hints_lib})
  set(mdsplus_LIBRARIES ${mdsplus_LIBRARY})
endif(mdsplus_USE_STATIC)

