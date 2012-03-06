from distutils.core import setup, Extension
 
module1 = Extension('mdsclient', 
        sources = ['mdsclientmodule.c'], 
        libraries = ['mdsclient'])
 
setup (name = 'mdsclient',
        version = '0.0',
        description = 'Mdsplus client',
        ext_modules = [module1])

