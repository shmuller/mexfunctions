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

    D = np.concatenate((D,b[:,None,:]),axis=1).copy()
    U = np.concatenate((U,np.zeros((J,1,M))),axis=1).copy()

    for j in range(J-1):
        U[j+1,-1] = D[j,-1]
        triblock.triblock(D[j],U[j+1],IPIV,L[j],D[j+1])

    b[J-1] = np.linalg.solve(D[J-1,:N].T,D[J-1,-1])
    for j in range(J-2,-1,-1):
        b[j] = U[j+1,-1] - np.dot(U[j+1,:N].T,b[j+1])


    """
    for j in range(J-1):
        triblock.triblock(D[j],U[j+1],IPIV,L[j],D[j+1],1)

    b[J-1] = np.linalg.solve(D[J-1].T,b[J-1])
    for j in range(J-2,-1,-1):
        b[j] = -np.dot(U[j+1].T,b[j+1])
    """

    return b

"""
if ~is_half
    b = reshape(b,[M,1,J]);
    D = [D,b];
    U = [U,zeros(M,1,J)];

    for j = 1:J-1
        U(:,N+1,j+1) = D(:,N+1,j);
        triblockmex(D(:,1:N,j),U(:,:,j+1),IPIV,L(:,:,j),D(:,:,j+1));
        %U(:,:,j+1) = D(:,1:N,j)\U(:,:,j+1);
        %D(:,1:N,j) = NaN;
        %D(:,:,j+1) = D(:,:,j+1) - L(:,:,j)*U(:,:,j+1);
    end

    b = reshape(b,[M,J]);
    b(:,J) = D(:,1:N,J)\D(:,N+1,J);
    for j = J-1:-1:1
        b(:,j) = U(:,N+1,j+1)-U(:,1:N,j+1)*b(:,j+1);
    end

else
    D = permute(D,[2,1,3]);
    L = permute(L,[2,1,3]);
    L = L(:,M/2+1:M,:);

    for j = 1:J-1
        triblockmex(D(:,:,j),U(:,:,j+1),IPIV,L(:,:,j),D(:,:,j+1),is_half);
    end

    b(:,J) = D(:,:,J).'\b(:,J);
    for j = J-1:-1:1
        b(:,j) = -U(:,:,j+1)*b(:,j+1);
    end


end

"""

