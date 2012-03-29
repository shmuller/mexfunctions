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

I2__ = dblMaxwmex('vrel',vt,v0,ut,u0,int32([3,3]),'MM')

v3 = linspace(-5,5,1000);

tic
I3 = dblMaxwmex('ion',vt,v0,ut,u0,int32([2,3]),'MM',v3);
toc

fun = @(x) dblMaxwmex('one',vt,v0,ut,u0,int32([2,3]),'MM',x);

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
I13_ion = dblMaxwmex('ion',vn,z,ve,z,int32([1,3]),'MM',v1,v2);
toc

figure;
surf(1e-3*v2,1e-3*v1,1e14*I13_ion)

tic
I13_CX = dblMaxwmex('CX',vn,z,vi,z,int32([1,3]),'MM',v1,v2);
toc

figure;
surf(1e-3*v2,1e-3*v1,1e14*I13_CX)


v3 = linspace(-5*vi,5*vi,100);
u3 = linspace(-5*vi,5*vi,50);

tic
I22_ion = dblMaxwmex('ion',vn,z,ve,z,int32([2,2]),'OO',v3,u3);
toc

figure;
surf(1e-3*u3,1e-3*v3,I22_ion);

tic
I22_CX = dblMaxwmex('CX',vn,z,vi,z,int32([2,2]),'OO',v3,u3);
toc

figure;
surf(1e-3*u3,1e-3*v3,I22_CX);

N1 = 100;
N2 = 50;
N3 = N1;
N4 = N2;
%N3 = 50;
%N4 = 100;

v1 = linspace(-5*vi,5*vi,N1).';
v2 = linspace(-5*vi,5*vi,N2).';
u1 = linspace(-5*vi,5*vi,N3).';
u2 = linspace(-5*vi,5*vi,N4).';

%v2 = 0;
%u2 = 0;

tic
I = dblMaxwmex('vrel',vn,z,vi,z,int32([1,1]),'MM',v1,v2,u1,u2);
toc

dv1 = diff(v1)/2; dv1 = [0;dv1] + [dv1;0];
dv2 = diff(v2)/2; dv2 = [0;dv2] + [dv2;0];
du1 = diff(u1)/2; du1 = [0;du1] + [du1;0];
du2 = diff(u2)/2; du2 = [0;du2] + [du2;0];

I = squeeze(sum(repmat(dv1,[1,N2,N3,N4]).*I,1));
I = squeeze(sum(repmat(dv2,[1,N3,N4]).*I,1));
I = sum(repmat(du1,[1,N4]).*I,1);
I = sum(du2.*I(:))

I_ana = 4/sqrt(2*pi)*hypot(vn,vi)


