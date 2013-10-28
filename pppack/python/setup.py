def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)

    """
    regenerate python wrappers with:
    > f2py pppack.f90 -m pppack -h pppack.pyf
    > mv pppack.pyf python
    > cd python

    then build with: sources=['pppack.pyf'], f2py_options=['--no-wrap-functions']

    to switch to pregenerated pppackmodule.c:
    > cp build/src.macosx-10.7-intel-2.7/* .

    and build again with no f2py_options and sources=['pppackmodule.c', 'fortranobject.c']
    """
    config.add_extension('pppack', 
                         #sources = ['pppack.pyf'],
                         #f2py_options = ['--no-wrap-functions'],
                         sources = ['pppackmodule.c', 'fortranobject.c'],
                         libraries = ['pppack'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

