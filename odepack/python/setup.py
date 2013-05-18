from distutils.core import setup, Extension
import numpy as np

module1 = Extension('odepack._odepack', 
        sources = ['odepack/odepackmodule.c'], 
        libraries = ['odepack'])
 
py_modules = [
        'odepack.odepack',
        'odepack.odepack_ctypes',
        'odepack.odepack_test']

setup (name = 'odepack',
        version = '0.0',
        description = 'Package for ODE solvers',
        include_dirs = [np.get_include()],
        py_modules = py_modules,
        ext_modules = [module1])

