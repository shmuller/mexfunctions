def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)

    """
    regenerate python wrappers with:
    > f2py *.f -m dierckx -h dierckx.pyf
    > vim dierckx.pyf
      :%s/real/real*8/g
      :x
    > mv dierckx.pyf python
    > cd python

    then build with: sources=['dierckx.pyf'], f2py_options=['--no-wrap-functions']

    to switch to pregenerated dierckxmodule.c:
    > cp build/src.macosx-10.7-intel-2.7/* .

    and build again with no f2py_options and sources=['dierckxmodule.c', 'fortranobject.c']
    """
    config.add_extension('dierckx', 
                         #sources = ['dierckx.pyf'],
                         #f2py_options = ['--no-wrap-functions'],
                         sources = ['dierckxmodule.c', 'fortranobject.c'],
                         libraries = ['dierckx'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

