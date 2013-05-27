from distutils.core import setup, Extension
from Cython.Distutils import build_ext

include_dirs = []
ext_modules = []
py_modules = ['fitfun_ctypes']

try:
    import fitfun_cffi
    ext_modules.append(fitfun_cffi.ffi.verifier.get_extension())
    py_modules.append('fitfun_cffi')
except ImportError:
    pass

try:
    import numpy as np
    include_dirs.append(np.get_include())

    module_accel = Extension('accel',
                             sources = ['accelmodule.c'])

    module_fitfun = Extension('fitfun', 
                              sources = ['fitfunmodule.c'],
                              libraries = ['fitfun'])

    from fitfun_cython_autogen import AutoGen

    funs = ('e2', 'IV3', 'IVdbl')
    funs_a = ('IV4', 'IV5', 'IV6', 'IVdbl2')

    autogen = AutoGen(funs=funs, funs_a=funs_a)
    autogen.pxd_gen()
    autogen.pyx_gen()

    module_fitfun_cython = Extension('fitfun_cython',
                                     sources = ['fitfun_cython.pyx'],
                                     libraries = ['fitfun'])

    ext_modules.extend([module_accel, module_fitfun, module_fitfun_cython])

except ImportError:
    # most likely pypy failing to import numpy, so install ctypes version
    # under fitfun.py
    py_modules.append('fitfun')
 
setup(name = 'fitfun',
      version = '0.0',
      description = 'Package for fitting functions',
      cmdclass = {'build_ext': build_ext},
      py_modules = py_modules,
      include_dirs = include_dirs,
      ext_modules = ext_modules)

