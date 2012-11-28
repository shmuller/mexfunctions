from distutils.core import setup, Extension
import numpy as np

module1 = Extension('atomic', 
        sources = ['atomicmodule.c'], 
        libraries = ['atomic'])
 
setup (name = 'atomic',
        version = '0.0',
        description = 'Package for atomic properities',
	include_dirs = [np.get_include()],
        ext_modules = [module1])

