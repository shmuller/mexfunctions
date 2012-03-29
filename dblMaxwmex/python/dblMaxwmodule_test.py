import numpy as np

import dblMaxw

v0 = np.array([0.1,0.2,0.3])
u0 = np.array([0.2,0.4,0.8])
vt = 0.5
ut = 2.

IJ = np.array([3,3],'i')

VU = ()

res = dblMaxw.integrate("vrel",vt,v0,ut,u0,IJ,VU)

print "%.16e" % res


qe = 1.60217646e-19
kB = 1.38065030e-23
mp = 1.67262158e-27
me = 9.10938188e-31

Te = 1500
ve = np.sqrt(qe/me*Te)

mi = 2*mp
Ti = 1500
vi = np.sqrt(qe/mi*Ti)

mn = 2*mp
Tn = 50
vn = np.sqrt(qe/mn*Tn)


N1, N2 = 100, 50
N3, N4 = N1, N2
#N3, N4 = 50, 100

v1 = np.linspace(-5*vi,5*vi,N1)
v2 = np.linspace(-5*vi,5*vi,N2)
u1 = np.linspace(-5*vi,5*vi,N3)
u2 = np.linspace(-5*vi,5*vi,N4)

z = np.zeros(3);
IJ = np.array([1,1],'i')

I = dblMaxw.integrate("vrel",vn,z,vi,z,IJ,(v1,v2,u1,u2))


dv1 = np.diff(v1)/2.; dv1 = np.r_[0.,dv1] + np.r_[dv1,0.]
dv2 = np.diff(v2)/2.; dv2 = np.r_[0.,dv2] + np.r_[dv2,0.]
du1 = np.diff(u1)/2.; du1 = np.r_[0.,du1] + np.r_[du1,0.]
du2 = np.diff(u2)/2.; du2 = np.r_[0.,du2] + np.r_[du2,0.]

I_sum = np.sum(dv1[np.newaxis,np.newaxis,np.newaxis,:]*I,3)
I_sum = np.sum(dv2[np.newaxis,np.newaxis,:]*I_sum,2)
I_sum = np.sum(du1[np.newaxis,:]*I_sum,1)
I_sum = np.sum(du2*I_sum,0)

I_ana = 4./np.sqrt(2*np.pi)*np.hypot(vn,vi)

print "I_sum = %.16e" % I_sum
print "I_ana = %.16e" % I_ana

