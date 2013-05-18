__all__ = ['odesolve', 'test_odesolve']

try:
    from odepack import odesolve
except:
    from odepack_ctypes import odesolve

from odepack_test import test_odesolve
