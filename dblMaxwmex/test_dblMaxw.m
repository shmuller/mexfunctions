v0 = [0.1,0.2,0.3];
u0 = [0.2,0.4,0.8];
vt = 0.5;
ut = 2;

v0(1:2) = 0;
u0(1:2) = 0;

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

I2__ = dblMaxwmex(vt,v0,ut,u0,int32([3,3]))


v3 = linspace(-5,5,1000);
tic
I3 = dblMaxwmex(vt,v0,ut,u0,int32([2,3]),v3);
toc

fun = @(x) dblMaxwmex(vt,v0,ut,u0,int32([2,3]),x);

tic
nrm = quadmex(int32(1),int32(256),fun,[-5,5])
toc
