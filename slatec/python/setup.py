def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)

    """
    regenerate python wrappers with:
    > f2py -m slatec -h slatec.pyf `cat f2py_files.txt`

    then build with: sources=['slatec.pyf'], f2py_options=['--no-wrap-functions']

    to switch to pregenerated slatecmodule.c:
    > cp build/src.macosx-10.7-intel-2.7/slatecmodule.c .

    and build again with no f2py_options and sources=['slatecmodule.c', 'fortranobject.c']
    """
    config.add_extension('slatec', 
                         #sources = ['slatec.pyf'],
                         #f2py_options = ['--no-wrap-functions'],
                         sources = ['slatecmodule.c', 'fortranobject.c'],
                         libraries = ['slatec'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

