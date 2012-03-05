import numpy as np
import triblock

def solve(L,D,U,b,is_half=0):
    """
    b = triblock_solve(L,D,U,b,is_half)
       Solve block-tridiagonal system Ax=b in place. L, D and U are 3D arrays 
       of N block matrices, such that

            | D1 U2               |
       A =  | L1 D2 ...           | 
            |           L(N-1) DN | 

       S. H. Muller, 2012/03/05
    """

    J, N, M = D.shape

    IPIV = np.zeros(M,dtype='i')
    
    if not is_half:
        D = np.concatenate((D,b[:,None,:]),axis=1).copy()
        U = np.concatenate((U,np.zeros((J,1,M))),axis=1).copy()

        for j in range(J-1):
            U[j+1,-1] = D[j,-1]
            triblock.step(D[j],U[j+1],IPIV,L[j],D[j+1])

        b[J-1] = np.linalg.solve(D[J-1,:N].T,D[J-1,-1])
        for j in range(J-2,-1,-1):
            b[j] = U[j+1,-1] - np.dot(U[j+1,:N].T,b[j+1])

    else:
        D = D.swapaxes(1,2).copy()
        L = L.swapaxes(1,2)[:,M/2:].copy()

        for j in range(J-1):
            triblock.step(D[j],U[j+1],IPIV,L[j],D[j+1],1)

        b[J-1] = np.linalg.solve(D[J-1],b[J-1])
        for j in range(J-2,-1,-1):
            b[j] = -np.dot(U[j+1].T,b[j+1])

    return b


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

    xx = np.c_[x,x2.ravel(),x3.ravel()]

    print xx


