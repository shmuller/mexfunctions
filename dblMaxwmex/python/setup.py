from distutils.core import setup, Extension
import numpy as np

module1 = Extension('dblMaxw', 
        sources = ['dblMaxwmodule.c'],
        libraries = ['dblMaxw'])
 
setup(name = 'dblMaxw',
        version = '0.0',
        description = 'Package for double-Maxwellian integrals',
	include_dirs = [np.get_include()],
        ext_modules = [module1])

