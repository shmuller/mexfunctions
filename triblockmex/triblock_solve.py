"""
    Solve block-tridiagonal system Ax=b in place. L, D and U are 3D arrays 
    of N block matrices, such that

        | D1 U2               |
    A = | L1 D2 ...           | 
        |           L(N-1) DN | 

    S. H. Muller, 2012/03/05
"""

import numpy as np
import triblock
import h5py

class save_restore:
    """Save and restore blocks."""
    def __init__(self, fname="tmp.h5", fmt="%03d", **kwargs):
        self.fname = fname
        self.fmt = fmt

        def donothing(j): pass
        self.print_save = kwargs.pop("print_save", donothing)
        self.print_load = kwargs.pop("print_load", donothing)

    def open(self, mode="w"):
        self.h = h5py.File(self.fname, mode)
    def close(self):
        self.h.close()
    def save(self, j, x):
        self.print_save(j)
        self.h[self.fmt % j] = x
    def load(self, j):
        self.print_load(j)
        return self.h[self.fmt % j][:]


class block_generator:
    """Generate blocks."""
    def __init__(self, L, D, U):
        M = D.shape[2]
        self.D = D.swapaxes(1,2).copy()
        self.L = L.swapaxes(1,2)[:,M/2:].copy()
        self.U = U.copy()

    def getL(self,j):
        return self.L[j]
    def getD(self,j):
        return self.D[j]
    def getU(self,j):
        return self.U[j]


def diagonalize(S,G,siz):
    J, M = siz
    IPIV = np.zeros(M,'i')

    D1 = G.getD(0)
    for j in range(J-1):
        D2 = G.getD(j+1)
        U2 = G.getU(j+1)
        L1 = G.getL(j)
        triblock.step(D1,U2,IPIV,L1,D2,1)
        D1 = D2
        S.save(j,U2)
    
    S.save(J-1,D1)


def backsubst(S,b,siz):
    J, M = siz
    D1 = S.load(J-1)
    b[J-1] = np.linalg.solve(D1,b[J-1])
    for j in range(J-2,-1,-1):
        U2 = S.load(j)
        b[j] = -np.dot(U2.T,b[j+1])
    return b


def solve_otf(S,G,b):
    siz = b.shape

    try:
        S.open("w")
        diagonalize(S,G,siz)
        b = backsubst(S,b,siz)
    finally:
        S.close()

    return b



def solve_full(L,D,U,b):
    """Solve full block-tridiagonal system."""
    J, N, M = D.shape
    IPIV = np.zeros(M,'i')

    D = np.concatenate((D,b[:,None,:]),axis=1).copy()
    U = np.concatenate((U,np.zeros((J,1,M))),axis=1).copy()

    for j in range(J-1):
        U[j+1,-1] = D[j,-1]
        triblock.step(D[j],U[j+1],IPIV,L[j],D[j+1])

    b[J-1] = np.linalg.solve(D[J-1,:N].T,D[J-1,-1])
    for j in range(J-2,-1,-1):
        b[j] = U[j+1,-1] - np.dot(U[j+1,:N].T,b[j+1])
    return b


def solve_half(L,D,U,b):
    """Solve full block-tridiagonal system."""
    J, N, M = D.shape
    IPIV = np.zeros(M,'i')

    D = D.swapaxes(1,2).copy()
    L = L.swapaxes(1,2)[:,M/2:].copy()
    U = U.copy()

    for j in range(J-1):
        triblock.step(D[j],U[j+1],IPIV,L[j],D[j+1],1)

    b[J-1] = np.linalg.solve(D[J-1],b[J-1])
    for j in range(J-2,-1,-1):
        b[j] = -np.dot(U[j+1].T,b[j+1])
    return b


def solve(L,D,U,b,is_half=0):
    """Solve block-tridiagonal system."""
    if not is_half:
        return solve_full(L,D,U,b)
    else:
        return solve_half(L,D,U,b)


if __name__ == "__main__":
    n, m = 2, 3

    D = np.random.rand(m,n,n)*5
    U = np.random.rand(m,n,n)*2
    L = np.random.rand(m,n,n)*3
    #U = np.zeros((m,n,n))
    #L = np.zeros((m,n,n))
    b = np.arange(1.,n*m+1.).reshape(m,n)

    U[:,:,n/2:] = 0.
    L[:,:,:n/2] = 0.
    b[:m-1,:] = 0.

    S = save_restore("tmp.h5", "%03d")
    G = block_generator(L,D,U)

    A = np.zeros((n*m,n*m))

    j = np.arange(n,dtype='i')
    for i in range(m):
        A[np.ix_(i*n+j,i*n+j)] = D[i]

    for i in range(1,m):
        A[np.ix_(i*n+j,(i-1)*n+j)] = U[i]

    for i in range(m-1):
        A[np.ix_(i*n+j,(i+1)*n+j)] = L[i]

    x = np.linalg.solve(A.T,b.ravel())

    x2 = solve(L,D,U,b.copy())

    x3 = solve(L,D,U,b.copy(),1)

    x4 = solve_otf(S,G,b.copy())

    xx = np.c_[x,x2.ravel(),x3.ravel(),x4.ravel()]

    print xx


