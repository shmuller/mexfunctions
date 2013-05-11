def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)

    config.add_extension('dierckx', 
                         sources = ['dierckx.pyf'], 
                         libraries = ['dierckx'], 
                         extra_f77_compile_args = ['-fdefault-real-8'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

