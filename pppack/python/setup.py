from distutils_f2py import ConfigurationF2Py

if __name__ == "__main__":
    #ConfigurationF2Py('pppack').build()
    #ConfigurationF2Py('dpppack', r8=True).build()
    #ConfigurationF2Py('pppack', from_pyf=False).build()
    ConfigurationF2Py('dpppack', from_pyf=False).build()

