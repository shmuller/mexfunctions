from distutils_f2py import ConfigurationF2Py

class ConfigurationSlatec(ConfigurationF2Py):
    def __init__(self, name, *args, **kw):
        ConfigurationF2Py.__init__(self, name, *args, **kw)
    
        self.add_extension('_' + name,
                           sources = ['_' + name + 'module.c'],
                           libraries = [name])


if __name__ == "__main__":
    #ConfigurationSlatec('slatec').build()
    #ConfigurationSlatec('dslatec', r8=True).build()
    ConfigurationSlatec('slatec', from_pyf=False).build()
    #ConfigurationSlatec('dslatec', from_pyf=False).build()

