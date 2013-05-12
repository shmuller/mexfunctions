def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)

    config.add_extension('dierckx', 
                         #sources = ['dierckx.pyf'],
                         #f2py_options = ['--no-wrap-functions'],
                         sources = ['dierckxmodule.c', 'fortranobject.c'],
                         libraries = ['dierckx'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

