from distutils.core import setup, Extension
 
module1 = Extension('mediansmooth', 
          sources = ['mediansmoothmodule.c'],
          libraries = ['mediansmooth'])
 
module2 = Extension('mediansmooth2', 
          sources = ['mediansmoothmodule.c'],
          libraries = ['mediansmooth2'])

setup(name = 'mediansmooth',
      version = '0.0',
      description = 'Package for median smoothing',
      ext_modules = [module1])

setup(name = 'mediansmooth2',
      version = '0.0',
      description = 'Package for median smoothing - pqueue2 backend',
      ext_modules = [module2])

