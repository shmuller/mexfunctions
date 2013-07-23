from distutils.core import setup, Extension

include_dirs = []
try:
    import numpy as np
except ImportError:
    pass
else:
    include_dirs.append(np.get_include())

module1 = Extension('mdsclient._mdsclient', 
                    sources=['mdsclient/mdsclientmodule.c'], 
                    libraries=['mdsclient'])

py_modules = ['mdsclient.__init__']

setup(name='mdsclient',
      version='0.0',
      description='Mdsplus client',
      include_dirs=include_dirs,
      py_modules=py_modules,
      ext_modules=[module1])

