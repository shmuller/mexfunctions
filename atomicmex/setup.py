from distutils.core import setup, Extension
 
module1 = Extension('atomic', 
        sources = ['atomicmodule.c', 'atomic.c'], 
        libraries = ['common'])
 
setup (name = 'atomic',
        version = '0.0',
        description = 'Package for atomic properities',
        ext_modules = [module1])

