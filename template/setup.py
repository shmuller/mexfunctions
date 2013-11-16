from distutils.core import setup, Extension
import numpy as np

setup(name='template',
      include_dirs=[np.get_include()],
      ext_modules=[Extension('template', sources=['templatemodule.c'])])

