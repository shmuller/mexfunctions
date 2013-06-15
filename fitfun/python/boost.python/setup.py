from distutils.core import setup, Extension
import numpy as np

extension_name = 'fitfun_boost'
extension_version = '0.0'

include_dirs = [np.get_include()]
library_dirs = []
libraries = ['boost_python-mt']
source_files = ['fitfun_boost.cpp']

# create the extension and add it to the python distribution
setup(name=extension_name, 
      version=extension_version, 
      ext_modules=[Extension(extension_name, 
                             source_files, 
                             include_dirs=include_dirs, 
                             library_dirs=library_dirs, 
                             libraries=libraries)])



extension_name = 'fitfun_boost_numpy'
extension_version = '0.0'

include_dirs = []
library_dirs = []
libraries = ['boost_python-mt', 'boost_numpy']
source_files = ['fitfun_boost_numpy.cpp']

# create the extension and add it to the python distribution
setup(name=extension_name, 
      version=extension_version, 
      ext_modules=[Extension(extension_name, 
                             source_files, 
                             include_dirs=include_dirs, 
                             library_dirs=library_dirs, 
                             libraries=libraries)])

