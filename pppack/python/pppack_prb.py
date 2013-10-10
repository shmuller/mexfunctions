import numpy as np
import pppack
from matplotlib.pyplot import plot, show

# test14
n = 12
ntitan = 49
k = 5
lenscr = (n-k)*(2*k+3)+5*k+3

a = np.zeros(n)

gtitan = np.zeros(ntitan)
gtau = np.zeros(ntitan)
i = 0
iflag = 0
ipick = np.array((1, 5, 11, 21, 27, 29, 31, 33, 35, 40, 45, 49), np.int32)
scrtch = np.zeros(lenscr)
t = np.zeros(n+k)
tau = np.zeros(n)
x = np.zeros(ntitan)

lx = pppack.titand(x, gtitan)

tau = x[ipick - 1]
gtau = gtitan[ipick - 1]

pppack.splopt(tau, k, scrtch, t, iflag)

pppack.splint(tau, gtau, t, k, scrtch[:(2*k-1)*n], a, iflag)

gtitan2 = np.zeros(ntitan)

for i in xrange(lx):
    gtitan2[i] = pppack.bvalue(t, a, k, x[i], 0)

#plot(x, gtitan, tau, gtau, 'r*', x, gtitan2)
#show()

#test17
k = 4
npoint = 201
brk = np.zeros(122)
coef = np.zeros((k, 22), order='F')
scrtch = np.zeros((n, 6), order='F')

l, k, iflag = pppack.tautsp(tau, gtau, 0., scrtch, brk, coef)

plott = np.zeros(npoint)
plotf = np.zeros(npoint)

step = (tau[-1] - tau[0]) / (npoint - 1)
for i in xrange(npoint):
    plott[i] = tau[0] + step * i
    plotf[i] = pppack.ppvalu(brk[:l+1], coef[:,:l], plott[i], 0)

gamma = 2.5

l, k, iflag = pppack.tautsp(tau, gtau, gamma, scrtch, brk, coef)

plotts = np.zeros(npoint)
for i in xrange(npoint):
    plotts[i] = pppack.ppvalu(brk[:l+1], coef[:,:l], plott[i], 0)

plot(x, gtitan, tau, gtau, 'r*', plott, np.c_[plotf, plotts])
show()


