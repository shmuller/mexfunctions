import numpy as np
import pppack
from matplotlib.pyplot import plot, show

class Bsplines:
    def __init__(self):
        self.t = np.array((0., 0., 0., 1., 1., 3., 4., 6., 6., 6.))
        
    def test04(self):
        k = 3
        t = self.t
        n = t.size - k
        
        npoint = 31*10

        xl = t[k-1]
        dx = (t[n] - xl) / (npoint - 1)
        x = xl + np.arange(npoint) * dx

        values = np.zeros((npoint, n+k))
        for i in xrange(npoint):
            left, mflag = pppack.interv(t[:n+1], x[i])
            pppack.bsplvb(t[:left+k], 1, x[i], left, values[i,left-k:left])
         
        plot(x, values[:,k-1:n])


bsplines = Bsplines()
bsplines.test04()


class Titan:
    def __init__(self):
        ntitan = 49        
        self.x = np.zeros(ntitan)
        self.gtitan = np.zeros(ntitan)
        
        lx = pppack.titand(self.x, self.gtitan)
        
        ipick = np.array((1, 5, 11, 21, 27, 29, 31, 33, 35, 40, 45, 49))

        self.tau = self.x[ipick - 1]
        self.gtau = self.gtitan[ipick - 1]

    def test14(self):
        n = 12
        k = 5
        lenscr = (n-k)*(2*k+3)+5*k+3
        scrtch = np.zeros(lenscr)
        t = np.zeros(n+k)

        iflag = pppack.splopt(self.tau, k, scrtch, t)

        a = np.zeros(n)

        iflag = pppack.splint(self.tau, self.gtau, t, k, scrtch[:(2*k-1)*n], a)

        x = self.x
        ntitan = x.size
        gtitan2 = np.zeros(ntitan)

        for i in xrange(ntitan):
            gtitan2[i] = pppack.bvalue(t, a, k, x[i], 0)

        plot(x, self.gtitan, self.tau, self.gtau, 'r*', x, gtitan2)

    def test17(self):
        k = 4
        npoint = 201
        brk = np.zeros(122)
        coef = np.zeros((k, 22), order='F')
        scrtch = np.zeros((n, 6), order='F')

        l, k, iflag = pppack.tautsp(self.tau, self.gtau, 0., scrtch, brk, coef)

        plott = np.zeros(npoint)
        plotf = np.zeros(npoint)

        tau0 = self.tau[0]
        tau1 = self.tau[-1]
        step = (tau1 - tau0) / (npoint - 1)
        for i in xrange(npoint):
            plott[i] = tau0 + step * i
            plotf[i] = pppack.ppvalu(brk[:l+1], coef[:,:l], plott[i], 0)

        gamma = 2.5

        l, k, iflag = pppack.tautsp(tau, gtau, gamma, scrtch, brk, coef)

        plotts = np.zeros(npoint)
        for i in xrange(npoint):
            plotts[i] = pppack.ppvalu(brk[:l+1], coef[:,:l], plott[i], 0)

        plot(x, gtitan, tau, gtau, 'r*', plott, np.c_[plotf, plotts])


titan = Titan()

#titan.test14()
#titan.test17()
show()

