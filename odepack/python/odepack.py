import numpy as np

import _odepack
from odepack_test import test_odesolve

def odesolve(f, y, t, ng=0, g=None, ibbox=None, bbox=None):
    n, neq = y.shape

    work = np.zeros(neq), np.zeros(1), np.zeros(neq), np.zeros(ng)

    if ng == 0:
        return _odepack.odesolve(f, y, t, work)
    else:
        return _odepack.odesolve(f, y, t, work, g, ibbox, bbox)


if __name__ == "__main__":
    test_odesolve(odesolve)

