function b = triblock_solve(L,D,U,b,is_half)
%b = triblock_solve(L,D,U,b)
%   Solve block-tridiagonal system Ax=b in place. L, D and U are 3D arrays 
%   of N block matrices, such that
%
%        | D1 U2               |
%   A =  | L1 D2 ...           | 
%        |           L(N-1) DN | 
%
%   S. H. Muller, 2011/12/12

if nargin < 5
    is_half = 0;
end

[M,N,J] = size(D);

IPIV = zeros(M,1,'int32');

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





%b(:,1) = D(:,:,1)\b(:,1);
%D(:,:,1) = D(:,:,1)\U(:,:,2);
%for j = 2:J-1
%    A = D(:,:,j)-L(:,:,j-1)*D(:,:,j-1);
%    b(:,j) = A\(b(:,j)-L(:,:,j-1)*b(:,j-1));
%    D(:,:,j) = A\U(:,:,j+1);
%end
%
%A = D(:,:,J)-L(:,:,J-1)*D(:,:,J-1);
%b(:,J) = A\(b(:,J)-L(:,:,J-1)*b(:,J-1));
%
%for j = J-1:-1:1
%    b(:,j) = b(:,j)-D(:,:,j)*b(:,j+1);
%end
