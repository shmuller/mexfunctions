from distutils.core import setup, Extension
import numpy as np

module1 = Extension('fieldline._fieldline', 
        sources = ['fieldline/fieldlinemodule.c'], 
        libraries = ['fieldline'])
 
py_modules = ['fieldline.fieldline']

setup (name = 'fieldline',
        version = '0.0',
        description = 'Package for field line integration',
        include_dirs = [np.get_include()],
        py_modules = py_modules,
        ext_modules = [module1])

