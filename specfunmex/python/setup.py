from distutils.core import setup, Extension
import numpy as np

module1 = Extension('specfun', 
        sources = ['specfunmodule.c'], 
        libraries = ['specfun', 'common'])
 
setup (name = 'specfun',
        version = '0.0',
        description = 'Package for special functions',
	include_dirs = [np.get_include()],
        ext_modules = [module1])

