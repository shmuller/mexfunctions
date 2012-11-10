from distutils.core import setup, Extension
import numpy as np

module1 = Extension('triblock', 
        sources = ['triblockmodule.c'],
        libraries = ['_triblock'])
 
setup(name = 'triblock',
        version = '0.0',
        description = 'Package for tri-block-diagonal matrices',
	include_dirs = [np.get_include()],
        ext_modules = [module1],
        py_modules = ['triblock_solve'])

