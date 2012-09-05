from distutils.core import setup, Extension
 
module1 = Extension('atomic', 
        sources = ['atomicmodule.c'], 
        libraries = ['atomic'])
 
setup (name = 'atomic',
        version = '0.0',
        description = 'Package for atomic properities',
        ext_modules = [module1])

