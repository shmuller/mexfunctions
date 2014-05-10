import numpy as np
find, cat, cont = np.flatnonzero, np.concatenate, np.ascontiguousarray
diff, zeros, ones, arange = np.diff, np.zeros, np.ones, np.arange
sign = np.sign

import operator

def Len(x):
    try:
        return len(x)
    except TypeError:
        return 1

from matplotlib.pyplot import plot, figure, show
from mpl_toolkits.mplot3d import Axes3D

def surf(x, y, Z, ax=None, cmap='jet'):
    if ax is None:
        fig = figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x[:,None], y[None,:], Z, rstride=1, cstride=1, cmap=cmap)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    return ax

def searchsorted(meshsites, sites):
    # need stable sorting algorithm to work properly
    index = cat((meshsites, sites)).argsort(kind='mergesort')
    return find(index > len(meshsites)-1) - np.arange(1, len(sites)+1)

#def searchsorted(meshsites, sites):
#    return meshsites.searchsorted(sites, 'right') - 1

def brk2knt(brks, mults):
    s = mults.sum()
    li = brks.size
    mm = np.zeros(s, 'i')
    mm[mults[:li-1].cumsum()] = 1
    t = brks[mm.cumsum()]
    return t

def augknt(knots, k, mults=np.array(1), with_addl=False):
    dk = np.diff(knots)
    if np.any(dk < 0):
        dk = np.diff(np.sort(knots))

    j = np.flatnonzero(dk)
    addl = k - j[0]
    interior = np.arange(j[0]+1, j[-1]+1)

    mults = cont(mults)
    if mults.size != interior.size:
        mults = mults[0].repeat(interior.size)

    augknot = brk2knt(knots[np.r_[0, interior, -1]], np.r_[k, mults, k])
    if not with_addl:
        return augknot
    else:
        return augknot, addl

def aveknt(t, k):
    n = t.size - k
    tstar = np.zeros(n)
    for j in xrange(n):
        tstar[j] = t[j+1:j+k].mean()
    return tstar

def aptknt_old(x, k=4):
    return augknt(cat((x[:1], aveknt(x, k), x[-1:])), k)

def aptknt(x, k=4):
    if x.dtype == np.float64:
        import dslatec as slatec
    else:
        import slatec
    t = np.zeros(x.size + k, x.dtype)
    slatec.aptknt(x, k, t)
    return t

def optknt(x, k=4):
    if k < 3:
        return aptknt(x, k)
    if x.dtype == np.float64:
        import dpppack as pppack
    else:
        import pppack
    n = x.size
    t = np.zeros(n+k)
    scrtch = np.zeros((n-k)*(2*k+3) + 5*k + 3)
    iflag = pppack.splopt(x, k, scrtch, t)
    return t

def chbpnt(t, k, tol=0.001, itermax=10):
    n = len(t) - k
    tau = aveknt(t, k)
    difftau = diff(tau)
    b = (arange(n - 1, -1, -1) % 2) * 2. - 1.
    sp = SplinePGS(tau, b, t=t)

    r = arange(1, n - 1)
    s = zeros(n - 2, np.int_)
    for _ in xrange(itermax):
        sp.plot()
        Dsp = sp.deriv()
        DDsp = Dsp.deriv()
        intau = tau[1:n-1]
        Dsptau = Dsp(intau).ravel()
        DDsptau = DDsp(intau).ravel()
        inb = b[1:n-1]
        sign(inb * Dsptau, s)
        dt = tau[r + s] - intau

        sign(inb * DDsptau, s)
        i = find(s >= 0)
        DDsptau[i] = (-2) * (2 * inb[i] / dt[i] + Dsptau[i]) / dt[i]
        dtau = -Dsptau / DDsptau
        tau[1:n-1] += dtau

        difftauold = difftau
        while True:
            difftau = diff(tau)
            if all(difftau > 0.1 * difftauold):
                break
            dtau /= 2.
            tau[1:n-1] -= dtau

        extremes = abs(sp(tau).ravel())
        sp = SplinePGS(tau, b, t=t)

        minext = min(extremes)
        maxext = max(extremes)
        if (maxext - minext) <= tol * minext:
            return tau, sp

    raise Exception("Failed to reach tolerance %f in %d iterations" % (tol, itermax))


def get_left(t, k, n, x):
    x = cont(x)
    left = t.searchsorted(x, 'right')
    left[left == 0] = k
    left[left == n+k] = n
    return left


class PP:
    def __init__(self, *args, **kw):
        self.ppmak(*args, **kw)

    def ppmak(self, b, a):
        self.b, self.a = b, a
        self.l = len(b) - 1
        self.p = a.shape[0]
        self.k = a.shape[1]
        self.d = a.shape[2] / self.l

    def ppbrk(self):
        return self.b, self.a, self.p, self.l, self.k, self.d

    def zeros(self):
        b, a, p, l, k, d = self.ppbrk()
        anew = np.zeros((p, 1, d*l))
        return self.__class__(b, anew)

    def _deriv(self, dorder):
        b, a, p, l, k, d = self.ppbrk()
        knew = k - dorder
        fact = np.ones(knew)
        expo = np.arange(k-1, 0, -1)
        for j in xrange(dorder):
            fact[:] *= expo[j:j+knew]
        anew = a[:,:knew] * fact[None,:,None]
        return self.__class__(b, anew)

    def deriv(self, dorder=1):
        if dorder < self.k:
            return self._deriv(dorder)
        else:
            return self.zeros()

    def ppual(self, x):
        if any(diff(x) < 0):
            tosort = True
            ix = np.argsort(x)
            xs = x[ix]
        else:
            xs = x.copy()

        b, a, p, l, k, d = self.ppbrk()

        index = self.get_index(x)
        xs[:] -= b[index]

        xs = xs[None,:,None]
        a = a.reshape((p, k, l, d))
        v = a[:,0,index].copy()
        for i in xrange(1, k):
            v[:] = xs * v + a[:,i,index]
        return v

    def get_index(self, sites):
        ind = searchsorted(self.b[:-1], sites)
        ind[ind < 0] = 0
        return ind

    def to_pp_pgs(self):
        b, a, p, l, k, d = self.ppbrk()
        anew = cont(a[:,::-1].reshape((p, k, l, d)).transpose((0, 2, 1, 3)))
        return PPPGS(b, anew)


class Spline(object):
    def __init__(self, *args, **kw):
        self.spmak(*args, **kw)
    
    def __call__(self, x):
        return self.spval(x)

    @classmethod
    def from_knots_coefs(cls, t, c):
        self = cls.__new__(cls)
        self.spmak(t, c)
        return self

    def spmak(self, t, c):
        self.t, self.c = cont(t, c.dtype), c
        self.p, self.n, self.d = c.shape
        self.k = t.size - self.n

    def spbrk(self):
        return self.t, self.c, self.k, self.p, self.n, self.d

    def _op_factory(op):
        def apply(self, other):
            return self.from_knots_coefs(self.t, op(self.c, other))
        return apply

    __add__ = _op_factory(operator.add)
    __sub__ = _op_factory(operator.sub)
    __mul__ = _op_factory(operator.mul)
    __div__ = _op_factory(operator.div)

    def bbox(self):
        return self.t[0], self.t[-1]

    def tave(self):
        t, k = self.t, self.k
        return aveknt(t, k)

    def get_left(self, x):
        return get_left(self.t, self.k, self.n, x)

    def spval(self, x, der=0, y=None):
        t, c, k, p, n, d = self.spbrk()
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), c.dtype)
        index = self.get_left(x) - 1
        if k == 1:
            y[:,:] = c[:,index]
        else:
            tx, i = self._setup_tx_i(t, x, k, d, index, backwd=False)
            c = c.reshape((p, n*d))
            y = y.reshape((p, m*d))
            for j in xrange(p):
                yj = c[j,i]
                self._sprval(tx, yj, k)
                y[j] = yj[0]
        return y.reshape((p, m, d))

    def eval(self, x, fill=None):
        shape = x.shape
        x = x.ravel()
        y = self.spval(x).ravel()
        if fill is not None:
            y[(x < self.t[0]) | (self.t[-1] < x)] = fill
        return y.reshape(shape)

    def zeros(self):
        t, c, k, p, n, d = self.spbrk()
        t = t[k-1:n+1]
        c = zeros((p, n-k+1, d), c.dtype)
        return self.from_knots_coefs(t, c)

    def _deriv(self, dorder):
        t, c, k, p, n, d = self.spbrk()
        for j in xrange(k-1, k-dorder-1, -1):
            tt = t[j+1:j+n] - t[1:n]
            i = find(tt > 0)
            temp = diff(c, axis=1)
            c = temp[:,i] * (j / tt[i])[None,:,None]
            t = cat((t[i+1], t[n+1:n+j+2]))
            n = len(i)
        return self.from_knots_coefs(t, c)

    def deriv(self, dorder=1):
        if dorder < self.k:
            return self._deriv(dorder)
        else:
            return self.zeros()

    def to_pp(self, backwd=True):
        t, c, k, p, n, d = self.spbrk()

        inter = find(diff(t) > 0)
        b = cat((t[inter], t[inter[-1]+1:inter[-1]+2]))
        if k == 1:
            a = c[:,inter]
        else:
            tx, i = self._setup_tx_i(t, t[inter], k, d, inter, backwd=backwd)
            if p == 1:
                a = c.ravel()[i]
                self._sprpp(tx, a, k, backwd=backwd)
            else:
                c = c.reshape((p, n*d))
                a = np.zeros((p,) + i.shape)
                for j in xrange(p):
                    a[j] = c[j,i]
                    self._sprpp(tx, a[j], k, backwd=backwd)
        return PP(b, a.reshape(p, k, -1))

    def plot(self):
        from matplotlib.pyplot import plot
        bbox = self.bbox()
        x = np.linspace(bbox[0], bbox[-1], 200)
        y = self(x)
        plot(x, y[0])

    def plot_ctrlpoly(self):
        from matplotlib.pyplot import plot
        plot(self.tave(), self.c[0])

    @staticmethod
    def _setup_tx_i(t, x, k, d, inter, backwd=False):
        if backwd:
            o = np.arange(k-1, 1-k, -1)
            kd = np.arange(0, -k*d, -d)
        else:
            o = np.arange(2-k, k)
            kd = np.arange(d*(1-k), d, d)

        tx = t[inter[None,:]+o[:,None]] - x[None,:]

        i = (d*(inter+1)-1)[None,:,None] + np.arange(1-d, 1)[None,None,:] \
          + kd[:,None,None]

        tx = tx[:,:,None].repeat(d, axis=2).reshape((2*(k-1),-1))
        i = i.reshape((k, -1))
        return tx, i

    @staticmethod
    def _sprval(tx, b, k):
        for r in xrange(1, k):
            for i in xrange(k-r):
                tx0 = tx[i+r-1]
                tx1 = tx[i+k-1]
                b[i] = (tx1 * b[i] - tx0 * b[i+1]) / (tx1 - tx0)

    @staticmethod
    def _sprval_bw(tx, b, k):
        for r in xrange(1, k):
            for j in xrange(k-1, r-1, -1):
                tx0 = tx[j-1+k-r]
                tx1 = tx[j-1]
                b[j] = (tx1 * b[j] - tx0 * b[j-1]) / (tx1 - tx0)

    @staticmethod
    def _sprdif(tx, b, k):
        for r in xrange(1, k):
            factor = float(k-r) / r
            for i in xrange(k-1, r-1, -1):
                b[i] = (b[i] - b[i-1]) * factor / tx[i+k-1-r]

    @staticmethod
    def _sprdif_bw(tx, b, k):
        for r in xrange(1, k):
            factor = float(k-r) / r
            for j in xrange(k-r):
                b[j] = (b[j] - b[j+1]) * factor / tx[j-1+r]

    @staticmethod
    def _sprpp(tx, b, k, backwd=True):
        if backwd:
            Spline._sprval_bw(tx, b, k)
            Spline._sprdif_bw(tx, b, k)
        else:
            Spline._sprval(tx, b, k)
            Spline._sprdif(tx, b, k)
            b[:] = b[::-1]


import pppack, dpppack

class PPPGS(PP):
    def ppmak(self, b, a):
        self.b, self.a = cont(b, a.dtype), a
        self.p, self.l, self.k, self.d = a.shape
        if a.dtype == np.float64:
            self.pppack = dpppack
        else:
            self.pppack = pppack

    def zeros(self):
        b, a, p, l, k, d = self.ppbrk()
        anew = np.zeros((p, l, 1, d), a.dtype)
        return self.__class__(b, anew)

    def _deriv(self, dorder):
        b, a, p, l, k, d = self.ppbrk()
        knew = k - dorder
        fact = np.ones(knew)
        expo = np.arange(k-1, 0, -1)
        for j in xrange(dorder):
            fact[:] *= expo[j:j+knew]
        anew = a[:,:,k-knew:] * fact[None,None,::-1,None]
        return self.__class__(b, anew)
  
    def ppual(self, x, der=0, fast=True):
        b, a, p, l, k, d = self.ppbrk()
        dtype = a.dtype
        x = cont(x, dtype)
        m = x.size
        y = np.zeros((p, m, d), dtype)
        if der == 0 and fast:
            self.pppack.ppual(b, a.T, x, y.T)
        else:
            ppualder = self.pppack.ppualder
            for j in xrange(p):
                ppualder(b, a[j].T, x, der, y[j].T)
        return y


class PPPGS2(PPPGS):
    def ppual(self, x, der=0):
        # this uses the PGS normalization
        b, a, p, l, k, d = self.ppbrk()
        dtype = a.dtype
        x = cont(x, dtype)
        m = x.size
        y = np.zeros((p, m, d), dtype)
        a = cont(a.transpose((0, 3, 1, 2)))
        ppvalu = self.pppack.ppvalu
        for j in xrange(p):
            yj = y[j]
            aj = a[j]
            for i in xrange(m):
                xi = x[i]
                yji = yj[i]
                for dd in xrange(d):
                    yji[dd] = ppvalu(b, aj[dd].T, xi, der)
        return y


class SplinePGS(Spline):
    def __init__(self, x, y, k=4, t=None, c=None, getknt=aptknt):
        x, y = self._checkargs(x, y)
        dtype = y.dtype
        p, n, d = y.shape
        if t is None:
            t = getknt(x, k)
        q = np.zeros((2*k-1)*n, dtype)
        g = np.zeros(n, dtype)
        bcoef = np.zeros(n, dtype)
        if c is None:
            c = np.zeros((p, n, d), dtype)
        if dtype == np.float64:
            splint = dpppack.splint
        else:
            splint = pppack.splint
        for j in xrange(p):
            yj = y[j]
            cj = c[j]
            for dd in xrange(d):
                g[:] = yj[:,dd]
                iflag = splint(x, g, t, k, q, bcoef)
                cj[:,dd] = bcoef
        self.spmak(t, c)

    def _checkargs(self, x, y):
        x = cont(x, y.dtype)
        if y.ndim == 1:
            y = y[None,:,None]
        elif y.ndim == 2:
            y = y[None,:,:]
        elif y.ndim > 3:
            raise ValueError("y with more than 3 dimensions not supported")
        return x, y

    def spmak(self, t, c):
        self.t, self.c = cont(t, c.dtype), c
        self.p, self.n, self.d = c.shape
        self.k = t.size - self.n
        if c.dtype == np.float64:
            self.pppack = dpppack
        else:
            self.pppack = pppack
        self.w_kxk = np.zeros((self.k, self.k), c.dtype)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        if der == 0 and fast:
            self.pppack.spual(t, c.T, k, x, y.T)
        else:
            self.pppack.spualder(t, c.T, k, x, y.T, der)
        return y

    def evalB(self, x, der=0):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        left = self.get_left(x)
        m = x.size
        B = np.zeros((m, k), dtype)
        if der == 0:
            bsplvb = self.pppack.bsplvb
            for i in xrange(m):
                bsplvb(t[:1], 1, x[i], left[i], B[i])
        else:
            a = self.w_kxk
            dbiatx = np.zeros((der+1, k), dtype)
            bsplvd = self.pppack.bsplvd
            for i in xrange(m):
                bsplvd(t[:1], x[i], left[i], a.T, dbiatx.T)
                B[i] = dbiatx[der]
        return left, B

    def spval2(self, x, der=0, y=None):
        t, c, k, p, n, d = self.spbrk()
        x = cont(x)
        left, B = self.evalB(x, der)
        m = x.size
        B = B[:,None,:,None]
        if y is None:
            y = np.zeros((p, m, d), c.dtype)
        for i in xrange(m):
            y[:,i,:] = (c[:,left[i]-k:left[i],:] * B[i]).sum(axis=1)
        return y

    spval3 = spval2

    def _deriv(self, dorder):
        t, c, k, p, n, d = self.spbrk()
        cnew = c.copy()
        self.pppack.spder(t, cnew.T, k, dorder)
        if p == 1:
            cnew.resize((1, n-dorder, d))
        else:
            cnew = cnew[:,:n-dorder].copy()
        t = t[dorder:n+k-dorder]
        return self.from_knots_coefs(t, cnew)

    def to_pp(self):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        l = len(np.unique(t)) - 1

        b = np.zeros(l+1, dtype)
        scrtch = np.zeros((k, k, p, d), dtype)
        a = np.zeros((p, l, k, d), dtype)

        self.pppack.bsplppd(t, c.T, scrtch.T, b, a.T)
        return PPPGS(b, a)

    def to_pp2(self):
        # this uses the PGS normalization
        t, c, k, p, n, d = self.spbrk()
        l = n+1-k

        b = np.zeros(l+1)
        scrtch = np.zeros((d, k, k))
        a = np.zeros((p, l, k, d))

        for j in xrange(p):
            cj = c[j].T.copy()
            dpppack.bspp2d(t, cj.T, scrtch.T, b, a[j].T)

        l = len(np.unique(t)) - 1
        if l < n+1-k:
            b.resize(l+1)
            if p == 1:
                a.resize((1, l, k, d))
            else:
                a = a[:,:l].copy()
        return PPPGS2(b, a)


import slatec, dslatec

class PPSLA2(PPPGS2):
    def ppmak(self, b, a):
        self.b, self.a = cont(b, a.dtype), a
        self.p, self.l, self.k, self.d = a.shape
        if a.dtype == np.float64:
            self.slatec = dslatec
        else:
            self.slatec = slatec

    def ppual(self, x, der=0):
        # this uses the PGS normalization
        b, a, p, l, k, d = self.ppbrk()
        dtype = c.dtype 
        x = cont(x, dtype)
        m = x.size
        y = np.zeros((p, m, d), dtype)
        if der >= k:
            return y
        a = cont(a.transpose((0, 3, 1, 2)))
        inppv = self.get_index(x) + 1
        ppval = self.slatec.ppval
        for j in xrange(p):
            yj = y[j]
            aj = a[j]
            for i in xrange(m):
                xi = x[i]
                yji = yj[i]
                inppvi = inppv[i]
                for dd in xrange(d):
                    yji[dd] = ppval(aj[dd].T, b, l, k, der, xi, inppvi)
        return y


class SplineSLA(SplinePGS):
    def spmak(self, t, c):
        SplinePGS.spmak(self, t, c)
        if c.dtype == np.float64:
            self.slatec = dslatec
        else:
            self.slatec = slatec
        self.ileft = np.ones(1, np.int32)
        self.w_3k = np.zeros(3*self.k, c.dtype)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        if der >= k:
            return y
        inbv = self.get_left(x).astype('i')
        work = self.w_3k
        ind = find((t[0] <= x) & (x <= t[-1]))
        c = cont(c.transpose((0, 2, 1)))
        bvalu = self.slatec.bvalu
        for j in xrange(p):
            yj = y[j]
            cj = c[j]
            for i in ind:
                xi = x[i]
                yji = yj[i]
                inbvi = inbv[i]
                for dd in xrange(d):
                    yji[dd] = bvalu(t, cj[dd], n, k, der, xi, inbvi, work)
        return y
    
    def evalB(self, x, der=0):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype 
        x = cont(x, dtype)
        m = x.size
        left = self.get_left(x)
        B = np.zeros((m, k), dtype)
        if der >= k:
            return left, B
        ind = find((t[0] <= x) & (x <= t[-1]))
        if der == 0:
            bspvn = self.slatec.bspvn
            work = self.w_3k[:2]
            iwork = self.ileft
            for i in ind:
                bspvn(t, k, k, 1, x[i], left[i], B[i], work, iwork)
        else:
            bspvd = self.slatec.bspvd
            a = np.zeros((k, k), dtype)
            dbiatx = np.zeros((der+1, k), dtype)
            work = np.zeros((k+1)*(k+2)/2, dtype)
            for i in ind:
                bspvd(t, k, der+1, x[i], left[i], dbiatx.T, work)
                B[i] = dbiatx[der]
        return left, B

    def _deriv(self, dorder):
        # This is a memory-wasting approach to spline differentiation, since
        # all the temporary derivatives are kept in 'ad'.
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        nderiv = dorder+1
        cnew = np.zeros((p, n-dorder, d), dtype)
        ad = np.zeros((2*n-nderiv+1)*nderiv/2, dtype)
        c = cont(c.transpose((0, 2, 1)))
        bspdr = self.slatec.bspdr
        for j in xrange(p):
            yj = y[j]
            cj = c[j]
            cjnew = cnew[j]
            for dd in xrange(d):
                bspdr(t, cj[dd], n, k, nderiv, ad)
                cjnew[:,dd] = ad[dorder-n:]
        t = t[dorder:n+k-dorder]
        return self.from_knots_coefs(t, cnew)

    def spval3(self, x, der=0, y=None):
        # 'dbspdr' does nothing else than differentiating the spline 'der' times, 
        # storing the coefficients one after each other in the linear array 'ad', 
        # where 'dbsped' retrieves them from the correct location. Then 'dbsped' 
        # evaluates the spline using the approach of spval2().
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype 
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        if der >= k:
            return y
        nderiv = der+1
        inev = self.ileft
        work = self.w_3k
        ad = np.zeros((2*n-nderiv+1)*nderiv/2, dtype)
        c = cont(c.transpose((0, 2, 1)))
        ind = find((t[0] <= x) & (x <= t[-1]))
        bspdr = self.slatec.bspdr
        bsped = self.slatec.bsped
        for j in xrange(p):
            yj = y[j]
            cj = c[j]
            for dd in xrange(d):
                bspdr(t, cj[dd], n, k, nderiv, ad)
                for i in ind:
                    bsped(t, ad, n, k, nderiv, x[i], inev, yj[i,dd:dd+1], work)
        return y

    def to_pp2(self):
        # this uses the PGS normalization
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        l = len(np.unique(t)) - 1
        b = np.zeros(l+1, dtype)
        scrtch = np.zeros((d, k, k), dtype)
        a = np.zeros((p, d, l, k), dtype)
        work = np.zeros(k*(n+3), dtype)

        c = cont(c.transpose((0, 2, 1)))
        bsppp = self.slatec.bsppp
        for j in xrange(p):
            cj = c[j]
            aj = a[j]
            for i in xrange(d):
                bsppp(t, cj[i], n, k, aj[i].T, b, l+1, work)
        a = cont(a.transpose((0, 2, 3, 1)))
        return PPSLA2(b, a)


class SplineSLA1(SplineSLA):
    def dbval(self, *args):
        return self.slatec.dbval1(*args)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros(m, dtype)
        inbv = self.ileft
        work = self.w_3k
        self.dbval(t, c.ravel(), k, der, x, inbv, work, y)
        return y.reshape((1, m, 1))


class SplineSLAI(SplineSLA1):
    def dbval(self, *args):
        return self.slatec.dbvali(*args)

import _slatec, _dslatec

class SplineSLAIC(SplineSLAI):
    def dbval(self, t, c, *args):
        if c.dtype == np.float64:
            return _dslatec.dbvali(t, c, *args)
        else:
            return _slatec.dbvali(t, c, *args)


class SplineSLA2(SplineSLA):
    def dbual(self, *args):
        return self.slatec.dbualu(*args)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        inbv = self.ileft
        work = np.zeros(k*(k+2), dtype)
        self.dbual(t, c.T, k, der, x, inbv, work, y.T)
        return y


class SplineSLA3(SplineSLA):
    def dbual(self, *args):
        return self.slatec.dbualu2(*args)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        inbv = self.ileft
        work = np.zeros(k*(k+1), dtype)
        work2 = np.zeros((k, d), dtype)
        self.dbual(t, c.T, k, der, x, inbv, work, work2.T, y.T)
        return y


class SplineSLA4(SplineSLA):
    def dbual(self, *args):
        return self.slatec.dbual(*args)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        inbv = self.ileft
        work = self.w_3k
        self.dbual(t, c.T, k, der, x, inbv, work, y.T)
        return y

class SplineSLA5(SplineSLA4):
    def dbual(self, *args):
        return self.slatec.dbual2(*args)


class SplineSLA6(SplineSLA):
    def dbual(self, *args):
        return self.slatec.dbual3(*args)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        inbv = self.ileft
        work = self.w_3k
        work2 = np.zeros((k, p, d), dtype)
        self.dbual(t, c.T, k, der, x, inbv, work, work2.T, y.T)
        return y


class SplineSLA7(SplineSLA4):
    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((1, m, p*d), dtype)
        inbv = self.ileft
        work = self.w_3k
        c = cont(c.transpose((1, 0, 2)).reshape((1, n, p*d)))
        self.dbual(t, c.T, k, der, x, inbv, work, y.T)
        return cont(y.reshape((m, p, d)).transpose((1, 0, 2)))


class SplineSLA8(SplineSLA):
    def dbual(self, *args):
        return self.slatec.dbual4(*args)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((m, p*d), dtype)
        inbv = self.ileft
        work = self.w_3k
        c = cont(c.transpose((1, 0, 2)).reshape((n, p*d)))
        self.dbual(t, c.T, k, der, x, inbv, work, y.T)
        return cont(y.reshape((m, p, d)).transpose((1, 0, 2)))


class SplineSLA9(SplineSLA):
    def dbual(self, *args):
        return self.slatec.dbualnd(*args)

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        if y is None:
            y = np.zeros((p, m, d), dtype)
        inbv = self.ileft
        work = self.w_3k

        n = np.array([n], 'i')
        k = np.array([k], 'i')
        s = np.ones(1, 'i')
        der = np.array([der], 'i')
        self.dbual(t, c.ravel(), n, k, s, der, x[:,None].T, inbv, work, y.ravel())
        return y


import dierckx, ddierckx

class SplineDie(SplinePGS):
    def spmak(self, t, c):
        SplinePGS.spmak(self, t, c)
        if c.dtype == np.float64:
            self.dierckx = ddierckx
        else:
            self.dierckx = dierckx

    def spval(self, x, der=0, y=None, fast=True):
        t, c, k, p, n, d = self.spbrk()
        dtype = c.dtype
        x = cont(x, dtype)
        m = x.size
        pd = p*d
        if y is None:
            y = np.zeros(m*pd, dtype)
        c = c.transpose((0, 2, 1)).ravel()
        ier = 0
        self.dierckx.splevv(t, c, k-1, x, y, pd, ier)
        return cont(y.reshape((m, p, d)).transpose((1, 0, 2)))


class SplineND(object):
    SplineClass = SplineSLA4
    def __init__(self, x, y, k=4):
        SplineClass = self.SplineClass
        n = np.array(y.shape, np.int32)
        c = y.copy()
        nd = c.ndim
        t = []
        k = np.array(k, np.int32)
        if k.size == 1:
            k = k.repeat(nd)
        for d in xrange(nd):
            c = c.reshape((n[:d].prod(), n[d], n[d+1:].prod()))
            sp = SplineClass(x[d], c, k[d], c=c)
            t.append(sp.t)
            c = sp.c
        self.spmak(t, c.reshape(n))

    def __call__(self, x, der=0):
        return self.spval(x, der=der)
        
    def get_left(self, x):
        return np.array(map(get_left, self.t, self.k, self.n, x))

    @classmethod
    def from_knots_coefs(cls, t, c):
        self = cls.__new__(cls)
        self.spmak(t, c)
        return self

    def spmak(self, t, c):
        self.t, self.c = t, c
        self.t = [cont(t, c.dtype) for t in self.t]
        self.n = np.array(c.shape, np.int32)
        self.k = np.array([t.size for t in self.t], np.int32) - self.n
        if c.dtype == np.float64:
            self.slatec = dslatec
        else:
            self.slatec = slatec
        self.s = np.zeros(c.ndim, np.int32)
        self.inbv = np.ones(c.ndim, np.int32)
        self.work = np.zeros(self.k.sum(), c.dtype)

    def reduce1(self, x):
        t, k = self.t, self.k
        nd = len(t)
        i = self.get_left(x)
        t = [t[d][i[d]-k[d]:i[d]+k[d]] for d in xrange(nd)]
        s = [slice(i[d]-k[d], i[d]) for d in xrange(nd)]
        return self.from_knots_coefs(t, cont(self.c[s]))

    def _spval_checkx(self, x, dtype):
        if isinstance(x, np.ndarray):
            return np.atleast_2d(cont(x, dtype))
        else:
            nd = len(x)
            m = max(map(Len, x))
            X = np.zeros((m, nd), dtype)
            for d in xrange(nd):
                X[:,d] = x[d]
            return X

    def spval(self, x, der=0, y=None):
        t, c, n, k = cat(self.t), self.c, self.n, self.k
        dtype = self.c.dtype
        x = self._spval_checkx(x, dtype)
        m, nd = x.shape
        der = np.array(der, np.int32)
        if der.size == 1:
            der = der.repeat(nd)
        if y is None:
            y = np.zeros(m, dtype)
        self.slatec.dbualnd(t, c.ravel(), n, k, self.s, der, x.T, self.inbv, self.work, y)
        return y

    def spval_dims(self, x, dims=None, der=0):
        t, c, n, SplineClass = self.t, self.c, self.n.copy(), self.SplineClass
        nd, dtype = c.ndim, c.dtype
        if not hasattr(x, "__iter__"):
            x = [x]
        x = [cont(xi, dtype) for xi in x]

        m = np.array(map(len, x))
        if dims is None:
            dims = np.arange(len(x))
        # process dimensions in order of max reduction of points
        perm = (m-n[dims]).argsort()
        der = np.array(der, np.int32)
        if der.size == 1:
            der = der.repeat(len(x))
        for i in perm:
            d = dims[i]
            c = c.reshape((n[:d].prod(), n[d], n[d+1:].prod()))
            sp = SplineClass.from_knots_coefs(t[d], c)
            c = sp.spval(x[i], der[i])
            n[d] = m[i]
        c = c.reshape(n)
        if nd == len(dims):
            return c
        else:
            t = [t[i] for i in xrange(nd) if i not in dims]
            return self.from_knots_coefs(t, c.squeeze())

    def spval_grid(self, x, der=0, y=None, fast=True):
        t, c, n, k = cat(self.t), self.c, self.n, self.k
        nd, dtype = c.ndim, c.dtype
        m = np.array(map(Len, x), np.int32)
        X = np.zeros(m.sum(), dtype)
        j = 0
        for i in xrange(m.size):
            X[j:j+m[i]] = x[i]
            j += m[i]

        der = np.array(der, np.int32)
        if der.size == 1:
            der = der.repeat(nd)

        i = np.ones(m.sum(), np.int32)
        B = np.zeros(k.dot(m), dtype)
        if y is None:
            y = np.zeros(m, dtype)

        if fast and nd == 3:
            self.slatec.dbual3d(t, c.ravel(), n, k, der, X, m, i, B, y.ravel())
        else:
            s = np.zeros(4*nd, np.int32)
            self.slatec.dbualgd(t, c.ravel(), n, k, s, der, X, m, i, B, y.ravel())
        return y

    def plot(self):
        t = self.t
        nd = len(t)
        if nd == 1:
            t = t[0]
            x = np.linspace(t[0], t[-1], 200)
            y = self(x[:,None])
            plot(x, y)
        elif nd == 2:
            tx, ty = t
            x = np.linspace(tx[0], tx[-1], 20)
            y = np.linspace(ty[0], ty[-1], 20)
            z = self.spval_grid((x, y))
            surf(x, y, z)
        else:
            raise NotImplementedError("cannot plot more than 2 dimensions")


try:
    mgc = get_ipython().magic
except NameError:
    pass


test1 = test2 = test3 = test4 = test5 = bench = False
if __name__ == "__main__":
    test1 = True
    #bench = True

if test1:
    p, n, d, k, m, der = 1, 16, 1, 4, 101, 1
    #c = np.zeros((p, n, d))
    #for i in xrange(d): c[:,n-1-i,i] = 1.
    c = cont(np.random.rand(p, n, d), np.float32)

    #knots = np.sort(np.random.rand(n-k+2.))
    knots = np.arange(n-k+2.)
    #knots = np.array((0., 1.2, 2.1, 2.1, 4.3, 5.2, 6.1, 7.))
    t = augknt(knots, k)
    sp = Spline.from_knots_coefs(t, c)

    x = np.linspace(knots[0], knots[-1], m)

    dsp = sp.deriv(der)

    y = dsp.spval(x)
    
    pp = sp.to_pp()
    dpp = pp.deriv(der)
    y2 = dpp.ppual(x)

    sp_pgs = SplineSLA4.from_knots_coefs(t, c)

    y3 = sp_pgs.spval(x, der=der)
    y3b = sp_pgs.spval2(x, der=der)
    y3c = sp_pgs.spval3(x, der=der)

    dsp_pgs = sp_pgs.deriv(der)

    y4 = dsp_pgs.spval(x)
    y4b = dsp_pgs.spval(x, fast=False)
    y4c = dsp_pgs.spval2(x)

    pp_pgs = sp_pgs.to_pp()
    y5 = pp_pgs.ppual(x, der=der)

    dpp_pgs = pp_pgs.deriv(der)
    y6 = dpp_pgs.ppual(x)

    pp_pgs2 = sp_pgs.to_pp2()
    y7 = pp_pgs2.ppual(x, der=der)

    pp_pgsb = pp.to_pp_pgs()
    y8 = pp_pgsb.ppual(x, der=der)

    assert np.allclose(pp_pgs.a, pp_pgsb.a, atol=1e-6)

    s = 0
    from matplotlib.pyplot import plot, show
    plot(x, np.c_[y[s], y2[s], y3[s], y4[s], y5[s], y6[s], y7[s], y8[s]], '.-')
    show()

if bench:
    p, n, d, k, m, der = 1, 1000, 1, 4, 1000000, 0
    c = cont(np.random.rand(p, n, d), np.float32)

    knots = np.arange(n-k+2.)
    t = augknt(knots, k)
    x = np.linspace(knots[0], knots[-1], m)
    #np.random.shuffle(x)

    SplineClasses = (SplineDie, SplinePGS, SplineSLA1, SplineSLAI, SplineSLAIC,
            SplineSLA2, SplineSLA3, SplineSLA4, SplineSLA5, 
            SplineSLA6, SplineSLA7, SplineSLA8, SplineSLA9)
    for SplineClass in SplineClasses:
        sp = SplineClass.from_knots_coefs(t, c)
        print SplineClass.__name__
        mgc(u'%timeit sp.spval(x, der=der)')

    sp_pgs = SplinePGS.from_knots_coefs(t, c)
    pp_pgs = sp_pgs.to_pp().deriv(der)
    print pp_pgs.__class__.__name__
    mgc(u'%timeit pp_pgs.ppual(x)')

if test2:
    def f(x, y):
        return np.cos(x[:,None]) * np.cos(y[None,:])

    kx = ky = 4
    x0, y0 = np.arange(5.), np.arange(6.)
    Z0 = cont(f(x0, y0), np.float32)
    
    sp = SplineND((x0, y0), Z0, k=4)
    tx, ty = sp.t
    coefs = sp.c

    Z0_check = sp.spval_grid((x0, y0))
    print 'Z0 check:', (Z0 - Z0_check).ptp()

    x = np.linspace(0., 4., 10)
    y = np.linspace(0., 5., 12)
    Z = f(x, y)

    der = (0, 0)
    Z1 = sp.spval_grid((x, y), der=der)

    from pytokamak.utils import splines
    tck = (cont(ty, 'd'), cont(tx, 'd'), cont(coefs.T.ravel(), 'd'))
    sp2 = splines.Spline2D(data=dict(tck=tck, degrees=(3, 3)))
    Z2 = sp2.deriv(der[1], der[0])(y, x).T

    print 'Z1 check:', (Z1 - Z2).ptp()

    ax = surf(x, y, Z1, cmap='jet')
    #ax = surf(x, y, Z, ax=ax, cmap='hot')
    show()

    pos = np.array([[1.2, 3.2], [2.2, 3.2], [3.2, 3.2]])
    Z3 = sp(pos)
    print Z3

    Z4 = sp.spval_grid((np.array([1.2,2.2,3.2]), [3.2]))
    print Z4.squeeze()

if test3:
    nx, ny = 20, 25
    kx, ky = 3, 3
    tx = augknt(np.arange(nx-kx+2.), kx)
    ty = augknt(np.arange(ny-ky+2.), ky)

    c = cont(np.random.rand(nx, ny), np.float32)

    sp = SplineND.from_knots_coefs((tx, ty), c)

    x = np.linspace(tx[0], tx[-1], 140)
    y = np.linspace(ty[0], ty[-1], 100)

    Z = sp.spval_grid((x, y))

    from pytokamak.utils import splines
    tck = (cont(ty, 'd'), cont(tx, 'd'), cont(c.T.ravel(), 'd'))
    sp2 = splines.Spline2D(data=dict(tck=tck, degrees=(ky-1, kx-1)))
    Z2 = sp2(y, x).T

    print (Z - Z2).ptp()

    from scipy.ndimage import map_coordinates

    coordinates = np.meshgrid(x+0.5, y+0.5)

    Z3 = map_coordinates(c, coordinates, order=2, prefilter=False).T

    i = y.size / 2

    from matplotlib.pyplot import plot, show
    plot(x, Z[:,i], x, Z2[:,i], x, Z3[:,i])
    show()

if test4:
    #"""
    from pytokamak.tokamak import overview
    AUG = overview.AUGOverview(29733, eqi_dig='EQH')
    R, z, psi_n = AUG.eqi.R, AUG.eqi.z, AUG.eqi.psi_n
    # do not convert to double
    t = psi_n.t
    psi_n = cont(psi_n.amp(psi_n._x), np.float32)

    mgc('%time sp = SplineND((t, z, R), psi_n, k=(2, 4, 4))')
    #"""
    pos = AUG.XPR.pos.t_gt(1.).compressed()
    #tzR = np.c_[pos.t[:,None], pos.x[:,::-1]]
    tzR = (pos.t, pos.x[:,1], pos.x[:,0])

    print "Evaluating at %d trajectory positions..." % pos.t.size
    mgc('%time y = sp(tzR, der=(0,0,0))')
    plot(pos.t, y)
    show()
    
    n = (t.size, z.size, R.size)
    m = np.prod(n)
    print "Evaluating at %d grid positions..." % m
    mgc('%time psi_n1 = sp.spval_grid((t, z, R))')
    maxerr = (psi_n1 - psi_n).ptp()
    assert maxerr < 1e-5, maxerr

    """
    TZR = np.zeros(n + (3,))
    TZR[:,:,:,0] = t[:,None,None]
    TZR[:,:,:,1] = z[None,:,None]
    TZR[:,:,:,2] = R[None,None,:]
    
    print "Evaluating at %d grid positions, using line algorithm..." % m
    mgc('%time psi_n1 = sp(TZR.reshape((-1, 3))).reshape(n)')
    maxerr = (psi_n1 - psi_n).ptp()
    assert maxerr < 1e-14, maxerr
    """

if test5:
    k = 4
    #knots = np.arange(8.)
    knots = np.array((0., 1.2, 2.1, 2.1, 4.3, 6.0, 6.1, 7.))
    
    #c = np.zeros(knots.size - 2 + k)
    #c[4] = 1
    c = cont(np.random.rand(knots.size - 2 + k), np.float32)

    sp = SplineND.from_knots_coefs([augknt(knots, k)], c)

    m = 101
    x = np.linspace(knots[0], knots[-1], m)

    #y = sp.spval1([[2.5]])
    #print y
    der = 1
    sp2 = SplinePGS.from_knots_coefs(augknt(knots, k), c[None,:,None]).deriv(der)

    y = sp(x[:,None], der=der)
    plot(x, y, x, sp2(x)[0,:,0], '.')
    show()

