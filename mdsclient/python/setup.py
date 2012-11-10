from distutils.core import setup, Extension
import numpy as np

module1 = Extension('mdsclient', 
        sources = ['mdsclientmodule.c'], 
        libraries = ['mdsclient'])
 
setup (name = 'mdsclient',
        version = '0.0',
        description = 'Mdsplus client',
	include_dirs = [np.get_include()],
        ext_modules = [module1])

