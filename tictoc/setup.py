from distutils.core import setup, Extension

setup(name='tictoc',
      ext_modules=[Extension('tictoc', sources=['tictocmodule.c'])])

