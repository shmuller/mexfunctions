import os, shutil
import numpy as np
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

class ConfigurationF2Py(Configuration):
    f2py_cmd_fmt = "f2py -m {name} -h {name}.pyf " \
                   "`cat f2py_files.txt` --overwrite-signature"

    def __init__(self, name, from_pyf=True, r8=False, **kw):
        Configuration.__init__(self, **kw)
        self._name = name

        if from_pyf:
            sources = [name + '.pyf']
        else:
            sources = [name + 'module.c', 'fortranobject.c']

        self.add_extension(name, 
                           sources = sources,
                           f2py_options = ['--no-wrap-functions'],
                           libraries = [name])
        if from_pyf:
            if r8:
                self.build = self._build_from_pyf_r8
            else:
                self.build = self._build_from_pyf
        else:
            self.build = self._build

    def rm_build(self):
        shutil.rmtree('build', True)

    def _build(self):
        self.rm_build()
        setup(**self.todict())

    def _build_from_pyf(self):
        os.system(self.f2py_cmd_fmt.format(name=self._name))
        self._build()
        os.system('cp build/src*/' + self._name + 'module.c .')

    def _build_from_pyf_r8(self):
        with open('.f2py_f2cmap', 'w') as f:
            f.write("{'real':{'':'double'}}")
        self._build_from_pyf()
        os.system('rm .f2py_f2cmap')

