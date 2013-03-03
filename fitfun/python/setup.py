from distutils.core import setup, Extension
import numpy as np

module1 = Extension('fitfun', 
                    sources = ['fitfunmodule.c'],
                    libraries = ['fitfun'])
 
setup(name = 'fitfun',
      version = '0.0',
      description = 'Package for fitting functions',
      include_dirs = [np.get_include()],
      ext_modules = [module1])

