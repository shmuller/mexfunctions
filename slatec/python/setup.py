import os, shutil
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(from_pyf=True, name='dslatec'):
    config = Configuration()
    """
    regenerate python wrappers with:
    > f2py -m slatec -h slatec.pyf `cat f2py_files.txt`

    then build with: sources=['slatec.pyf'], f2py_options=['--no-wrap-functions']

    to switch to pregenerated slatecmodule.c:
    > cp build/src.macosx-10.7-intel-2.7/slatecmodule.c .

    and build again with no f2py_options and sources=['slatecmodule.c', 'fortranobject.c']
    """
    if from_pyf:
        sources = [name + '.pyf']
    else:
        sources = [name + 'module.c', 'fortranobject.c'],

    config.add_extension(name, 
                         sources = sources,
                         f2py_options = ['--no-wrap-functions'],
                         libraries = [name])
    
    config.add_extension('_slatec',
                         sources = ['_slatecmodule.c'],
                         libraries = ['dslatec'])
    return config


def build_double_from_pyf():
    print "*** Building from *.pyf with system fortranobject.c ***"
    shutil.rmtree('build', True)
    with open('.f2py_f2cmap', 'w') as f:
        f.write("{'real':{'':'double'}}")
    os.system('f2py -m dslatec -h dslatec.pyf `cat f2py_files.txt` --overwrite-signature')
    setup(**configuration(from_pyf=True, name='dslatec').todict())
    os.system('cp build/src*/dslatecmodule.c .')
    os.system('rm .f2py_f2cmap')

def build_double_from_module():
    print "*** Building from *module.c with local fortranobject.c ***"
    shutil.rmtree('build', True)
    setup(**configuration(from_pyf=False, name='dslatec').todict())


def build_single_from_pyf():
    print "*** Building from *.pyf with system fortranobject.c ***"
    shutil.rmtree('build', True)
    os.system('f2py -m slatec -h slatec.pyf `cat f2py_files.txt` --overwrite-signature')
    setup(**configuration(from_pyf=True, name='slatec').todict())
    os.system('cp build/src*/slatecmodule.c .')

def build_single_from_module():
    print "*** Building from *module.c with local fortranobject.c ***"
    shutil.rmtree('build', True)
    setup(**configuration(from_pyf=False, name='slatec').todict())


if __name__ == "__main__":    
    #build_double_from_pyf()
    #build_double_from_module()
    #build_single_from_pyf()
    build_single_from_module()

