from distutils.core import setup, Extension
 
module1 = Extension('dblMaxw', 
        sources = ['dblMaxwmodule.c'],
        libraries = ['dblMaxw'])
 
setup(name = 'dblMaxw',
        version = '0.0',
        description = 'Package for double-Maxwellian integrals',
        ext_modules = [module1])

