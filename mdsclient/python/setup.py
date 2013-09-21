from distutils.core import setup, Extension
import numpy as np

module1 = Extension('mdsclient._mdsclient', 
                    sources=['mdsclient/mdsclientmodule.c'], 
                    libraries=['mdsclient'])

py_modules = ['mdsclient.__init__']

setup(name='mdsclient',
      version='0.0',
      description='Mdsplus client',
      include_dirs=[np.get_include()],
      py_modules=py_modules,
      ext_modules=[module1])

