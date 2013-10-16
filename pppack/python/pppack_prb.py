import math
import numpy as np
import numpy.ma as ma
import pppack
from matplotlib.pyplot import figure, plot, show

def ev(tau, c, x):
    n = tau.size
    npoint = x.size
    y = np.zeros(npoint)
    for i in xrange(npoint):
        y[i] = pppack.ppvalu(tau, c[:,:n-1], x[i], 0)
    return y


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


#bsplines = Bsplines()
#bsplines.test04()
#bsplines.test05()
#bsplines.test06()
#bsplines.test07()

class Knots:
    def __init__(self):
        #self.itermx = 0
        self.itermx = 2
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

    def errmax_decay(self, tau, c, n):
        step = 20
        l = tau.size - 1

        errmax = 0.
        for i in xrange(l):
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

        figure()
        X = np.linspace(-1., 1., 1001)
        plot(X, self.g(X))

        for n in xrange(self.nlow, self.nhigh + 1, 2):
            h = 1. / (n - 1)
            tau = 2. * (h * np.arange(n))**self.irate09 - 1.

            def interp():
                c[0,:n] = self.g(tau)
                pppack.cubspl(tau, c[:,:n], 0, 0)            
            
            interp()
            for iter in xrange(self.itermx):
                pppack.newnot(tau.copy(), c[:,:n-1], tau, scrtch[:,:n-1])
                interp()    
                
            plot(X, ev(tau, c, X))

            self.errmax_decay(tau, c, n)
            
    def test10(self):
        k = 4
        bcoef = np.zeros(22)

        figure()
        X = np.linspace(-1., 1., 1001)
        plot(X, self.g(X))

        for n in xrange(self.nlow, self.nhigh + 1, 2):
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

            plot(X, ev(brk[:l+1], c[:,:l], X))

            self.errmax_decay(brk[:l+1], c[:,:l], n)

    def test12(self):
        nmax = 20
        k = 4
        bcoef = np.zeros(nmax + 2)
        scrtch = np.zeros((2*k-1)*nmax)
        scrtch2 = scrtch[:k*k].reshape((k, k), order='F')
        scrtch3 = scrtch[:2*nmax].reshape((2, nmax), order='F')
        brk = np.zeros(nmax)
        c = np.zeros((k, nmax), order='F')
        t = np.zeros(nmax + 6)
        tnew = np.zeros(nmax)
        tau = np.zeros(nmax)
        gtau = np.zeros(nmax)

        def interp(t, n):
            for i in xrange(n):
                tau[i] = t[i+1:i+4].sum() / 3.
            gtau[:n] = self.g(tau[:n])

            iflag = pppack.splint(tau[:n], gtau[:n], t[:n+k], k, scrtch[:(2*k-1)*n], bcoef[:n])
            l = pppack.bsplpp(t[:n+k], bcoef[:n], scrtch2, brk, c)
            return l

        figure()
        X = np.linspace(-1., 1., 1001)
        plot(X, self.g(X))

        for n in xrange(self.nlow, self.nhigh + 1, 2):
            if n == self.nlow:
                nmk = n - k
                h = 2. / (nmk + 1)
                t[:k] = -1
                t[k:n] = np.arange(1, nmk + 1) * h - 1.
                t[n:n+k] = 1.
            else:
                pppack.newnot(brk[:l+1], c[:,:l], tnew[:l+3], scrtch3[:,:l])
                l += 2
                t[4+l] = t[5+l] = 1.
                t[4:3+l] = tnew[1:l]
                        
            l = interp(t, n)
            for iter in xrange(self.itermx):
                pppack.newnot(brk[:l+1], c[:,:l], tnew[:l+1], scrtch3[:,:l])
                t[4:3+l] = tnew[1:l]
                l = interp(t, n)

            plot(X, ev(brk[:l+1], c[:,:l], X))

            self.errmax_decay(brk[:l+1], c[:,:l], n)


#knots = Knots()
#knots.test09()
#knots.test10()
#knots.test12()

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
        print "Optimal knots:", t[k:n]

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


#titan = Titan()
#titan.test14()
#titan.test17()

class L2Main:
    def __init__(self):
        ntmax = 200
        self.tau = np.zeros(ntmax)
        self.gtau = np.zeros(ntmax)
        self.weight = np.zeros(ntmax)
        self.brk = np.zeros(100)
        self.coef = np.zeros(2000)
        
        self.t = np.zeros(ntmax)
        self.scrtch = np.zeros(ntmax)

        self.setdat = pppack.setdatex4

    def ex2(self):
        icount = 0
        ntimes = 5
        addbrk = 2

        ntau, totalw, l, k = self.setdat(icount, self.tau, self.gtau, self.weight,
                                         self.brk, self.coef)
        tau = self.tau[:ntau]
        gtau = self.gtau[:ntau]

        t = self.t[:2*k-1+l]
        n = pppack.l2knts(self.brk[:l+1], k, t)
       
        def l2appr(t, k, n):
            q = np.zeros((k, n), order='F')
            bcoef = np.zeros(n)
            pppack.l2appr(t, q, self.scrtch[:n], bcoef, 
                          ntau, self.tau, self.gtau, self.weight)

            scrtch = q[:,:k]
            coef = self.coef[:k*n].reshape((k, n), order='F')
            l = pppack.bsplpp(t, bcoef, scrtch, self.brk, coef)
            return l, self.brk[:l+1], coef[:,:l]

        def l2err(l, k):
            prfun = 0
            ftau = np.zeros(ntau)
            error = np.zeros(ntau)
            pppack.l2err(prfun, ftau, error, 
                         self.tau, self.gtau, self.weight, totalw,
                         self.brk, self.coef, l, k)
            return ftau, error

        l, brk, coef = l2appr(t, k, n)
        ftau, error = l2err(l, k)

        lbegin = l
        for nt in xrange(ntimes):
            lnew = lbegin + nt * addbrk
            tnew = self.scrtch[:lnew+1]
            scrtch = self.t[:2*l].reshape((2, l), order='F')
            pppack.newnot(brk, coef, tnew, scrtch)

            t = self.t[:2*k-1+lnew]
            n = pppack.l2knts(tnew, k, t)

            l, brk, coef = l2appr(t, k, n)
            ftau, error = l2err(l, k)

        print brk

        x = np.linspace(tau[0], tau[-1], 1001)
        y = ev(brk, coef, x)

        figure()
        plot(tau, gtau, '-+')
        plot(tau, ftau)
        plot(x, y)
        plot(brk, np.zeros(brk.size), 'kd')

        figure()
        plot(tau, error, brk, np.zeros(brk.size), 'kd')


#l2main = L2Main()
#l2main.ex2()

class TensorProductSpline:
    def __init__(self):
        kx, nx = 3, 7
        taux = np.arange(1., nx + 1)
        
        ky, ny = 4, 6
        tauy = np.arange(1., ny + 1)
        
        work1 = np.zeros((ny, nx), order='F')
        work2 = np.zeros(max(nx, ny))
        work3 = np.zeros(max((2*kx-1)*nx, (2*ky-1)*ny))
        
        self.param = taux, tauy, kx, ky
        self.wrk = work1, work2, work3

    def g(self, x, y):
        return np.where(x < 3.5, 0., x-3.5)**2 + np.where(y < 3., 0., y-3.)**3

    def get_t(self, tau, k, mid=True):
        n = tau.size
        t = np.zeros(n+k)
        t[:k] = tau[0]
        t[n:] = tau[n-1]
        if mid:
            for i in xrange(k, n):
                t[i] = (tau[i-k+1] + tau[i-k+2]) / 2.
        else:
            for i in xrange(k, n):
                t[i] = tau[i-k+2]
        return t

    def spline(self, taux, tauy, bcoef, kx, ky):
        nx, ny = taux.size, tauy.size
        tx = self.get_t(taux, kx, mid=True)
        ty = self.get_t(tauy, ky, mid=False)
        work1, work2, work3 = self.wrk

        iflag = pppack.spli2d(taux, bcoef, tx, kx, work2[:nx], work3[:(2*kx-1)*nx], work1)
        iflag = pppack.spli2d(tauy, work1, ty, ky, work2[:ny], work3[:(2*ky-1)*ny], bcoef)
        return tx, ty, bcoef

    def ev(self, tx, ty, bcoef, x, y):
        nx, ny = bcoef.shape
        kx, ky = tx.size - nx, ty.size - ny
        mx, my = x.size, y.size
        res = np.zeros((my, mx))
        wrk = self.wrk[1]

        for j in xrange(my):
            lefty, mflag = pppack.interv(ty, y[j])
            for i in xrange(mx):
                for jj in xrange(ky):
                    wrk[jj] = pppack.bvalue(tx, bcoef[:,lefty-ky+jj], kx, x[i], 0)
                res[j,i] = pppack.bvalue(ty[lefty-ky:lefty+ky], wrk[:ky], ky, y[j], 0)
        return res

    def test19(self):
        taux, tauy, kx, ky = self.param
        
        bcoef = self.g(taux[None,:], tauy[:,None]).T
        
        tx, ty, bcoef = self.spline(taux, tauy, bcoef, kx, ky)

        fun = self.g(taux[None,:], tauy[:,None])

        res = self.ev(tx, ty, bcoef, taux, tauy)
        err = res - fun
        return fun, res, err


tps = TensorProductSpline()
fun, res, err = tps.test19()

show()

