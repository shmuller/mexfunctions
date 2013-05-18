from distutils.core import setup, Extension
import numpy as np

module1 = Extension('_odepack', 
        sources = ['odepackmodule.c'], 
        libraries = ['odepack'])
 
setup (name = 'odepack',
        version = '0.0',
        description = 'Package for ODE solvers',
        include_dirs = [np.get_include()],
        ext_modules = [module1])

