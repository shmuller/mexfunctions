n = 2;
m = 3;

D = rand(n,n,m)*5;
U = rand(n,n,m)*2;
L = rand(n,n,m)*3;
b = reshape(1:n*m,n,m);

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

