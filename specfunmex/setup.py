from distutils.core import setup, Extension
 
module1 = Extension('specfun', 
        sources = ['specfunmodule.c'], 
        libraries = ['common', 'specfun'])
 
setup (name = 'specfun',
        version = '0.0',
        description = 'Package for special functions',
        ext_modules = [module1])

