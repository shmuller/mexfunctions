from distutils.core import setup, Extension

include_dirs = []
try:
    import numpy as np
except ImportError:
    pass
else:
    include_dirs.append(np.get_include())

module1 = Extension('specfun', 
                    sources = ['specfunmodule.c'], 
                    libraries = ['specfun', 'common'])
 
setup(name = 'specfun',
      version = '0.0',
      description = 'Package for special functions',
      include_dirs = include_dirs,
      ext_modules = [module1])

