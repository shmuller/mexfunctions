m = 50;
phi = (0:m-1)*(2*pi/m);
cph = cos(phi); sph = sin(phi);

n1 = 1;
x1 = 2*rand(n1,1);
y1 = 2*rand(n1,1);
% x1 = [1,1,2,2,1].';
% y1 = [1,1,1,1,1].';
r1 = (1:n1).'./n1;

r1 = r1(end:-1:1);

[xy,n] = deal([]);
for j = 1:n1
    xy = [xy, [x1(j)+r1(j)*cph; y1(j)+r1(j)*sph]];
    n = [n; m];
end
c = [0,0,1; 0,1,0; 1,0,0; 0,0,0; 1,1,1];
%c = jet(n1);
a = 0.7*ones(n1,1);

gpc_tpatch(xy,n,c,a);

xlim([0,3]);
ylim([0,3]);
axis square;

