function b = triblock_solve(L,D,U,b)
%b = triblock_solve(L,D,U,b)
%   Solve block-tridiagonal system Ax=b in place. L, D and U are 3D arrays 
%   of N block matrices, such that
%
%        | D1 U2               |
%   A =  | L1 D2 ...           | 
%        |           L(N-1) DN | 
%
%   S. H. Muller, 2011/12/12

J = size(D,3);

b(:,1) = D(:,:,1)\b(:,1);
D(:,:,1) = D(:,:,1)\U(:,:,2);
for j = 2:J-1
    A = D(:,:,j)-L(:,:,j-1)*D(:,:,j-1);
    b(:,j) = A\(b(:,j)-L(:,:,j-1)*b(:,j-1));
    D(:,:,j) = A\U(:,:,j+1);
end

A = D(:,:,J)-L(:,:,J-1)*D(:,:,J-1);
b(:,J) = A\(b(:,J)-L(:,:,J-1)*b(:,J-1));

for j = J-1:-1:1
    b(:,j) = b(:,j)-D(:,:,j)*b(:,j+1);
end