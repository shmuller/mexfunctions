from distutils.core import setup, Extension

include_dirs = []
try:
    import numpy as np
except ImportError:
    pass
else:
    include_dirs.append(np.get_include())

module1 = Extension('sort', 
          sources = ['sortmodule.c'],
          libraries = ['sort'])
 
setup(name = 'sort',
      version = '0.0',
      description = 'Comparison of sorting algorithms',
      include_dirs = include_dirs,
      ext_modules = [module1])
