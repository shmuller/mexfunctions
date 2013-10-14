import math
import numpy as np
import numpy.ma as ma
import pppack
from matplotlib.pyplot import figure, plot, show

class Bsplines:
    def __init__(self):
        #original data for test04
        #self.k = 3
        #self.t = np.array((0., 0., 0., 1., 1., 3., 4., 6., 6., 6.))
        
        #original data for test05 and test06
        self.k = 4
        self.t = np.array((0., 0., 0., 0., 1., 3., 4., 6., 6., 6., 6.))
        figure()

    def test04(self):
        t, k = self.t, self.k
        n = t.size - k
        
        npoint = 31*10
        x = np.linspace(t[k-1], t[n], npoint)
        values = np.zeros((npoint, n))
        for i in xrange(npoint):
            left, mflag = pppack.interv(t[:n+1], x[i])
            pppack.bsplvb(t[:left+k], 1, x[i], left, values[i,left-k:left])
         
        plot(x, values)
        assert np.allclose(values.sum(axis=1), 1.)

    def test05(self):
        t, k = self.t, self.k

        npoint = 40
        x = (np.arange(npoint) + 1) / 5. - 1.
        lefts = np.arange(k, 2*k)
        biatx = np.zeros(k)
        values = np.zeros((npoint, k))
        for i in xrange(npoint):
            for left in lefts:
                pppack.bsplvb(t[:left+k], 1, x[i], left, biatx)
                values[i,left-k] = biatx[2*k-1-left]

        masked_values = ma.masked_array(values, np.abs(values) > 1.)
        plot(x, masked_values, '--')

    def test06(self):
        t, k = self.t, self.k
        n = t.size - k
        bcoef = np.zeros(n)
        bcoef[(n-1)/2] = 1.
        
        brk = np.zeros(n+1)
        coef = np.zeros((k, n), order='F')
        scrtch = np.zeros((k, k), order='F')

        l = pppack.bsplpp(t, bcoef, scrtch, brk, coef)

        npoint = 40
        x = (np.arange(npoint) + 1) / 5. - 1.
        values = np.zeros(npoint)
        for i in xrange(npoint):
            values[i] = pppack.ppvalu(brk[:l+1], coef[:,:l], x[i], 0)
        plot(x, values, 'k', linewidth=2)

    def test07(self):
        k = self.k
        t = self.t[k-1:1-k]
        npoint = 40
        x = (np.arange(npoint) + 1) * 0.2 - 1
        bcoef = np.ones(1)
        values = np.zeros(npoint)
        for i in xrange(npoint):
            values[i] = pppack.bvalue(t, bcoef, k, x[i], 0)
        plot(x, values, 'c--', linewidth=2)


bsplines = Bsplines()
bsplines.test04()
bsplines.test05()
bsplines.test06()
bsplines.test07()

class Knots:
    def __init__(self):
        self.itermx = 0
        #self.itermx = 3
        self.nlow = 4
        self.nhigh = 20

        self.irate09 = 1
        self.irate10 = 8

        def g(x):
            return np.sqrt(x + 1.)

        def dg(x):
            return 0.5 / g(x)

        def ddg(x):
            return -0.5 * dg(x) / (x + 1.)

        self.g, self.dg, self.ddg = g, dg, ddg

    def eval(self, tau, c, x):
        n = tau.size
        npoint = x.size
        y = np.zeros(npoint)
        for i in xrange(npoint):
            y[i] = pppack.ppvalu(tau, c[:,:n-1], x[i], 0)
        return y

    def errmax_decay(self, tau, c):
        step = 20
        n = tau.size

        errmax = 0.
        for i in xrange(n-1):
            dx = (tau[i+1] - tau[i]) / step
            for j in xrange(step):
                h = (j + 1) * dx
                x = tau[i] + h
                pnatx = c[0,i] + h * (c[1,i] + h * (c[2,i] + h * c[3,i] / 3.) / 2.)
                
                errmax = max(errmax, np.abs(self.g(x) - pnatx))

        decay = 0.
        aloger = np.log(errmax)
        if (n > self.nlow):
            decay = (aloger - self.algerp) / np.log(float(n) / float(n-2))
        self.algerp = aloger

        print n, errmax, decay

    def test09(self):
        nmax = 20
        c = np.zeros((4, nmax), order='F')
        scrtch = np.zeros((2, nmax), order='F')

        X = np.linspace(-1., 1., 1001)

        for n in xrange(self.nlow, self.nhigh + 1, 2):
            figure()
            h = 1. / (n - 1)
            tau = 2. * (h * np.arange(n))**self.irate09 - 1.

            c[0,:n] = self.g(tau)
            pppack.cubspl(tau, c[:,:n], 0, 0)            
            plot(X, self.g(X), X, self.eval(tau, c, X))

            for iter in xrange(self.itermx):
                pppack.newnot(tau.copy(), c[:,:n-1], tau, scrtch[:,:n-1])
                c[0,:n] = self.g(tau)
                pppack.cubspl(tau, c[:,:n], 0, 0)
                plot(X, self.eval(tau, c, X))

            self.errmax_decay(tau, c)
            
    def test10(self):
        bcoef = np.zeros(22)

        X = np.linspace(-1., 1., 1001)
        figure()
        plot(X, self.g(X))

        for n in xrange(self.nlow, self.nhigh + 1, 2):
            k = 4
            m = n + 2
            h = 1. / (n - 1)
            t = 2. * (h * np.arange(1, n - 1))**self.irate10 - 1.
            t = np.r_[[-1.] * k, t, [1.] * k]
        
            bcoef[0] = 0.
            dtip2 = t[4] - t[0]
            taui = t[4]
            bcoef[1] = self.g(taui) - 2. * dtip2 * self.dg(taui) / 3. \
                     + dtip2**2 * self.ddg(taui) / 6.
            
            for i in xrange(2, m):
                taui = t[i+2]
                dtip1 = dtip2
                dtip2 = t[i+3] - t[i+2]

                bcoef[i] = self.g(taui) + (dtip2 - dtip1) * self.dg(taui) / 3. \
                         - dtip1 * dtip2 * self.ddg(taui) / 6.

            brk = np.zeros(m+1)
            c = np.zeros((k, m), order='F')
            scrtch = np.zeros((k, k), order='F')

            l = pppack.bsplpp(t, bcoef[:m], scrtch, brk, c)

            plot(X, self.eval(brk[:l+1], c[:,:l], X))

            self.errmax_decay(brk[:l+1], c[:,:l])


knots = Knots()
knots.test09()
knots.test10()

class Titan:
    def __init__(self):
        ntitan = 49        
        self.x = np.zeros(ntitan)
        self.gtitan = np.zeros(ntitan)
        
        lx = pppack.titand(self.x, self.gtitan)
        
        ipick = np.array((1, 5, 11, 21, 27, 29, 31, 33, 35, 40, 45, 49))

        self.tau = self.x[ipick - 1]
        self.gtau = self.gtitan[ipick - 1]
        figure()

    def test14(self):
        n = self.tau.size
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
        n = self.tau.size
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

        l, k, iflag = pppack.tautsp(self.tau, self.gtau, gamma, scrtch, brk, coef)

        plotts = np.zeros(npoint)
        for i in xrange(npoint):
            plotts[i] = pppack.ppvalu(brk[:l+1], coef[:,:l], plott[i], 0)

        plot(self.x, self.gtitan, self.tau, self.gtau, 'r*', plott, np.c_[plotf, plotts])


titan = Titan()

titan.test14()
titan.test17()
show()

