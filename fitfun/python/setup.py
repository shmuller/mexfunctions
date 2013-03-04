from distutils.core import setup, Extension

include_dirs = []
ext_modules = []
py_modules = ['fitfun_ctypes']

try:
    import numpy as np
    include_dirs.append(np.get_include())

    module1 = Extension('fitfun', 
                        sources = ['fitfunmodule.c'],
                        libraries = ['fitfun'])
    ext_modules.append(module1)
except ImportError:
    # most likely pypy failing to import numpy, so install ctypes version
    # under fitfun.py
    py_modules.append('fitfun')
 
setup(name = 'fitfun',
      version = '0.0',
      description = 'Package for fitting functions',
      py_modules = py_modules,
      include_dirs = include_dirs,
      ext_modules = ext_modules)

