from distutils.core import setup, Extension
 
module1 = Extension('mediansmooth', 
        sources = ['mediansmoothmodule.c'],
        libraries = ['mediansmooth'])
 
setup(name = 'mediansmooth',
        version = '0.0',
        description = 'Package for median smoothing',
        ext_modules = [module1])

