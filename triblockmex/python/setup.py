from distutils.core import setup, Extension

include_dirs = []
try:
    import numpy as np
except ImportError:
    pass
else:
    include_dirs.append(np.get_include())

module1 = Extension('triblock', 
                    sources = ['triblockmodule.c'],
                    libraries = ['_triblock'])
 
setup(name = 'triblock',
      version = '0.0',
      description = 'Package for tri-block-diagonal matrices',
      include_dirs = include_dirs,
      ext_modules = [module1],
      py_modules = ['triblock_solve'])

