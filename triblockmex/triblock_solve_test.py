import numpy as np
import triblock_solve 

reload(triblock_solve)

n, m = 2, 3

D = np.random.rand(m,n,n)*5
U = np.random.rand(m,n,n)*2
L = np.random.rand(m,n,n)*3
#U = np.zeros((m,n,n))
#L = np.zeros((m,n,n))
b = np.arange(1.,n*m+1.).reshape(m,n)

#U[:,:,n/2:] = 0.
#L[:,:,:n/2] = 0.
#b[:m-1,:] = 0.

A = np.zeros((n*m,n*m))

j = np.arange(n,dtype='i')
for i in range(m):
    A[np.ix_(i*n+j,i*n+j)] = D[i,:,:].T

for i in range(1,m):
    A[np.ix_((i-1)*n+j,i*n+j)] = U[i,:,:].T

for i in range(m-1):
    A[np.ix_((i+1)*n+j,i*n+j)] = L[i,:,:].T

x = np.linalg.solve(A,b.reshape(-1))

x2 = triblock_solve.solve(L,D,U,b.copy())

xx = np.c_[x,x2.reshape(-1)]

print xx

"""
U(n/2+1:n,:,:) = 0;
L(1:n/2,:,:) = 0;
b(:,1:m-1) = 0;

% use standard solver
A = zeros(n*m,n*m);

for i = 1:m
    A((i-1)*n+(1:n),(i-1)*n+(1:n)) = D(:,:,i);
end

for i = 2:m
    A((i-2)*n+(1:n),(i-1)*n+(1:n)) = U(:,:,i);
end

for i = 1:m-1
    A(i*n+(1:n),(i-1)*n+(1:n)) = L(:,:,i);
end

x = A\b(:);

% mex solver
x2 = triblock_solve(L,D,U,b);

% mex solver with is_half
x3 = triblock_solve(L,D,U,b,1);


[x(:),x2(:),x3(:)]

"""
