n = 2;
m = 3;

D = rand(n,n,m)*5;
U = rand(n,n,m)*2;
L = rand(n,n,m)*3;
b = reshape(1:n*m,n,m);

x = triblock_solve(L,D,U,b);

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

x2 = A\b(:);

[x(:),x2(:)]

