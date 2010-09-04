function gpc_tpatch_test

m = 20;
phi = (0:m-1)*(2*pi/m);
cph = cos(phi); sph = sin(phi);

n1 = 10;
x1 = rand(n1,1);
y1 = rand(n1,1);
r1 = rand(n1,1);

[xy,N] = deal([]);
for j = 1:n1
    xy = [xy, [x1(j)+r1(j)*cph; y1(j)+r1(j)*sph]];
    N = [N; m];
end
xyca1 = {xy, N, [0,0,1], 0.5};

n2 = 15;
x2 = rand(n2,1);
y2 = rand(n2,1);
r2 = rand(n2,1);

[xy,N] = deal([]);
for j = 1:n2
    xy = [xy, [x2(j)+r2(j)*cph; y2(j)+r2(j)*sph]];
    N = [N; m];
end
xyca2 = {xy, N, [0,1,0], 0.5};

gpc_tpatch({xyca1,xyca2});

xlim([-1,2]);
ylim([-1,2]);
axis square;