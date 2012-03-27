v0 = [0.1,0.2,0.3];
u0 = [0.2,0.4,0.8];
vt = 0.5;
ut = 2;

%v0(1:2) = 0;
%u0(1:2) = 0;

L = [-5,5];

I = quadfunmex(int32([1,2,3]),int32([128,128,128]),'Maxw', ...
    [-10,10],[-10,10],[-10,10],v0(1),v0(2),v0(3),vt)


I_ = quadfunmex(int32(1),int32(256),'Maxw_r',[0,10],v0(1),v0(2),v0(3),vt)

tic
I2 = quadfunmex(int32([1,2,3,8,9,10]),int32([16,16,16,16,16,16]),'dblMaxw', ...
    L*vt,L*vt,L*vt,v0(1),v0(2),v0(3),vt, ...
    L*ut,L*ut,L*ut,u0(1),u0(2),u0(3),ut)
toc

w0 = u0-v0;
wt = hypot(vt,ut);
I2_ = quadfunmex(int32(1),int32(256),'Maxw_r',[0,10],w0(1),w0(2),w0(3),wt)

I2__ = dblMaxwmex('vrel',vt,v0,ut,u0,int32([3,3]))

v3 = linspace(-5,5,1000);

tic
I3 = dblMaxwmex('ion',vt,v0,ut,u0,int32([2,3]),v3);
toc

fun = @(x) dblMaxwmex('one',vt,v0,ut,u0,int32([2,3]),x);

tic
nrm = quadmex(int32(1),int32(256),fun,[-5,5])
toc



qe = 1.60217646e-19;
kB = 1.38065030e-23;
mp = 1.67262158e-27;
me = 9.10938188e-31;

Te = 1500;
ve = sqrt(qe/me*Te);

mi = 2*mp;
Ti = 1500;
vi = sqrt(qe/mi*Ti);

mn = 2*mp;
Tn = 50;
vn = sqrt(qe/mn*Tn);


v1 = linspace(-5*vi,5*vi,100);
v2 = linspace(-5*vi,5*vi,50);

z = [0,0,0];

tic
I13_ion = dblMaxwmex('ion',vn,z,ve,z,int32([1,3]),v1,v2);
toc

figure;
surf(1e-3*v2,1e-3*v1,1e14*I13_ion)

tic
I13_CX = dblMaxwmex('CX',vn,z,vi,z,int32([1,3]),v1,v2);
toc

figure;
surf(1e-3*v2,1e-3*v1,1e14*I13_CX)


v3 = linspace(-5*vi,5*vi,100);
u3 = linspace(-5*vi,5*vi,50);

tic
I22_ion = dblMaxwmex('ion',vn,z,ve,z,int32([2,2]),v3,u3);
toc

figure;
surf(1e-3*u3,1e-3*v3,I22_ion);

tic
I22_CX = dblMaxwmex('CX',vn,z,vi,z,int32([2,2]),v3,u3);
toc

figure;
surf(1e-3*u3,1e-3*v3,I22_CX);

v1 = linspace(-5*vi,5*vi,100);
u1 = linspace(-5*vi,5*vi,100);
v2 = linspace(-5*vi,5*vi,50);
u2 = linspace(-5*vi,5*vi,50);

%v2 = 0;
%u2 = 0;

tic
I = dblMaxwmex('CX',vn,z,vi,z,int32([1,1]),v1,u1,v2,u2);
toc


