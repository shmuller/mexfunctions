def configuration(from_pyf=True):
    from numpy.distutils.misc_util import Configuration
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
        sources = ['slatec.pyf']
    else:
        sources = ['slatecmodule.c', 'fortranobject.c'],

    config.add_extension('slatec', 
                         sources = sources,
                         f2py_options = ['--no-wrap-functions'],
                         libraries = ['slatec'])
    return config

if __name__ == "__main__":
    import os
    os.system('rm -rf build')

    from numpy.distutils.core import setup
    os.system('f2py -m slatec -h slatec.pyf `cat f2py_files.txt` --overwrite-signature')
    
    print "*** Building from *.pyf with system fortranobject.c ***"
    setup(**configuration(from_pyf=True).todict())
     
    os.system('cp build/src*/slatecmodule.c .')
    print "*** Building from *module.c with local fortranobject.c ***"
    os.system('rm -rf build/lib*/* build/temp*/*')
    setup(**configuration(from_pyf=False).todict())
