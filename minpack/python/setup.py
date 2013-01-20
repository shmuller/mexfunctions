from distutils.core import setup, Extension
import numpy as np

module1 = Extension('minpack', 
        sources = ['minpackmodule.c'], 
        libraries = ['minpack'])
 
setup (name = 'minpack',
        version = '0.0',
        description = 'Package for minimization functions',
        include_dirs = [np.get_include()],
        ext_modules = [module1])

