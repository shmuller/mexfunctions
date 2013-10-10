import numpy as np
import pppack

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

lx = ntitan
pppack.titand(x, gtitan, lx)

tau = x[ipick - 1]
gtau = gtitan[ipick - 1]

pppack.splopt(tau, k, scrtch, t, iflag)

pppack.splint(tau, gtau, t, k, scrtch[:(2*k-1)*n], a, iflag)

gtitan_ = np.zeros(ntitan)

for i in xrange(lx):
    gtitan_[i] = pppack.bvalue(t, a, k, x[i], 0)



