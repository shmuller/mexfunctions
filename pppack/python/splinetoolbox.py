import numpy as np
diff, find, cat, zeros, ones = \
        np.diff, np.flatnonzero, np.concatenate, np.zeros, np.ones

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

    mults = np.atleast_1d(mults)
    if mults.size != interior.size:
        mults = mults[0].repeat(interior.size)

    augknot = brk2knt(knots[np.r_[0, interior, -1]], np.r_[k, mults, k])
    if not with_addl:
        return augknot
    else:
        return augknot, addl


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

    def deriv(self, dorder=1):
        b, a, p, l, k, d = self.ppbrk()
        knew = k - dorder
        if knew > 0:
            fact = np.ones(knew)
            expo = np.arange(k-1, 0, -1)
            for j in xrange(dorder):
                fact[:] *= expo[j:j+knew]
            anew = a[:,:knew] * fact[None,:,None]
        else:
            anew = np.zeros((p, 1, d*l))
        return self.__class__(b, anew)

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
        anew = np.ascontiguousarray(a[:,::-1].reshape((p, k, l, d)).transpose((0, 2, 1, 3)))
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
        self.t, self.c = t, c
        self.p, self.n, self.d = c.shape
        self.k = t.size - self.n

    def spbrk(self):
        return self.t, self.c, self.k, self.p, self.n, self.d

    def bbox(self):
        return self.t[0], self.t[-1]

    def tave(self):
        t, k, n = self.t, self.k, self.n
        tstar = np.zeros(n)
        for j in xrange(n):
            tstar[j] = t[j+1:j+k].mean()
        return tstar

    def get_left(self, x):
        t, k, n = self.t, self.k, self.n
        left = t.searchsorted(x, 'right')
        left[left == n+k] = n
        return left

    def spval(self, x):
        t, c, k, p, n, d = self.spbrk()
        m = x.size
        index = self.get_left(x) - 1
        if k == 1:
            y = c[:,index]
        else:
            tx, i = self._setup_tx_i(t, x, k, d, index, backwd=False)
            if p == 1:
                y = c.ravel()[i]
                self._sprval(tx, y, k)
                y.resize((m, d))
            else:
                c = c.reshape((p, n*d))
                y = np.zeros((p, m*d))
                for j in xrange(p):
                    yj = c[j,i]
                    self._sprval(tx, yj, k)
                    y[j] = yj[0]
        return y.reshape((p, m, d))

    def deriv(self, dorder=1):
        t, c, k, p, n, d = self.spbrk()

        knew = k - dorder;
        if knew <= 0:
            t = t[k-1:n+1]
            c = zeros((p, n-k+1, d))
        else:
            for j in xrange(k-1, knew-1, -1):
                tt = t[j+1:j+n] - t[1:n]
                i = find(tt > 0)
                temp = diff(c, axis=1)
                c = temp[:,i] * (j / tt[i])[None,:,None]
                t = cat((t[i+1], t[n+1:n+j+2]))
                n = len(i)
        return self.from_knots_coefs(t, c)

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


import pppack

class PPPGS(PP):
    def ppmak(self, b, a):
        self.b, self.a = b, a
        self.p, self.l, self.k, self.d = a.shape

    def deriv(self, dorder=1):
        b, a, p, l, k, d = self.ppbrk()
        knew = k - dorder
        if knew > 0:
            fact = np.ones(knew)
            expo = np.arange(k-1, 0, -1)
            for j in xrange(dorder):
                fact[:] *= expo[j:j+knew]
            anew = a[:,:,k-knew:] * fact[None,None,::-1,None]
        else:
            anew = np.zeros((p, l, 1, d))
        return self.__class__(b, anew)

    def ppual(self, x, der=0, fast=True):
        b, a, p, l, k, d = self.ppbrk()

        m = x.size
        y = np.zeros((p, m, d))
        if der == 0 and fast:
            pppack.ppual(b, a.T, x, y.T)
        else:
            for j in xrange(p):
                pppack.ppualder(b, a[j].T, x, der, y[j].T)
        return y


class PPPGS2(PPPGS):
    def ppual(self, x, der=0):
        # this uses the PGS normalization
        b, a, p, l, k, d = self.ppbrk()
        m = x.size
        y = np.zeros((p, m, d))
        a = np.ascontiguousarray(a.transpose((0, 3, 1, 2)))
        for j in xrange(p):
            yj = y[j]
            aj = a[j]
            for i in xrange(m):
                xi = x[i]
                yji = yj[i]
                for dd in xrange(d):
                    yji[dd] = pppack.ppvalu(b, aj[dd].T, xi, der)
        return y


class SplinePGS(Spline):
    def spval(self, x, der=0, fast=True):
        t, c, k, p, n, d = self.spbrk()
        m = x.size
        y = np.zeros((p, m, d))
        if der == 0 and fast:
            pppack.spual(t, c.T, k, x, y.T)
        else:
            pppack.spualder(t, c.T, k, x, y.T, der)
        return y

    def evalB(self, x, der=0):
        t, c, k, p, n, d = self.spbrk()
        left = self.get_left(x)
        m = x.size
        B = np.zeros((m, k))
        if der == 0:
            for i in xrange(m):
                pppack.bsplvb(t[:1], 1, x[i], left[i], B[i])
        else:
            a = np.zeros((k, k))
            dbiatx = np.zeros((der+1, k))
            for i in xrange(m):
                pppack.bsplvd(t[:1], x[i], left[i], a.T, dbiatx.T)
                B[i] = dbiatx[der]
        return left, B

    def spval2(self, x, der=0):
        left, B = self.evalB(x, der)
        c, k, m = self.c, self.k, x.size
        B = B[:,None,:,None]
        y = np.zeros((self.p, m, self.d))
        for i in xrange(m):
            y[:,i,:] = (c[:,left[i]-k:left[i],:] * B[i]).sum(axis=1)
        return y

    def deriv(self, dorder=1):
        t, c, k, p, n, d = self.spbrk()
        if k <= dorder:
            t = t[k-1:n+1]
            cnew = zeros((p, n-k+1, d))
        else:
            cnew = c.copy()
            pppack.spder(t, cnew.T, k, dorder)
            if p == 1:
                cnew.resize((1, n-dorder, d))
            else:
                cnew = cnew[:,:n-dorder].copy()
            t = t[dorder:n+k-dorder]
        return self.from_knots_coefs(t, cnew)

    def to_pp(self):
        t, c, k, p, n, d = self.spbrk()
        l = len(np.unique(t)) - 1

        b = np.zeros(l+1)
        scrtch = np.zeros((k, k, p, d))
        a = np.zeros((p, l, k, d))

        pppack.bsplppd(t, c.T, scrtch.T, b, a.T)
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
            pppack.bspp2d(t, cj.T, scrtch.T, b, a[j].T)

        l = len(np.unique(t)) - 1
        if l < n+1-k:
            b.resize(l+1)
            if p == 1:
                a.resize((1, l, k, d))
            else:
                a = a[:,:l].copy()
        return PPPGS2(b, a)


import slatec

class PPSLA2(PPPGS2):
    def ppual(self, x, der=0):
        # this uses the PGS normalization
        b, a, p, l, k, d = self.ppbrk()
        m = x.size
        y = np.zeros((p, m, d))
        a = np.ascontiguousarray(a.transpose((0, 3, 1, 2)))
        inppv = self.get_index(x) + 1
        for j in xrange(p):
            yj = y[j]
            aj = a[j]
            for i in xrange(m):
                xi = x[i]
                yji = yj[i]
                inppvi = inppv[i]
                for dd in xrange(d):
                    yji[dd] = slatec.dppval(aj[dd].T, b, l, k, der, xi, inppvi)
        return y


class SplineSLA(SplinePGS):
    def spval(self, x, der=0, fast=True):
        t, c, k, p, n, d = self.spbrk()
        m = x.size
        y = np.zeros((p, m, d))
        inbv = self.get_left(x)
        work = np.zeros(3*k)
        c = np.ascontiguousarray(c.transpose((0, 2, 1)))
        for j in xrange(p):
            yj = y[j]
            cj = c[j]
            for i in xrange(m):
                xi = x[i]
                yji = yj[i]
                inbvi = inbv[i]
                for dd in xrange(d):
                    yji[dd] = slatec.dbualu(t, cj[dd], n, k, der, xi, inbvi, work)
        return y

    def evalB(self, x, der=0):
        t, c, k, p, n, d = self.spbrk()
        left = self.get_left(x)
        m = x.size
        B = np.zeros((m, k))
        work = np.zeros(2*k)
        iwork = np.zeros(1, 'i')
        if der == 0:
            for i in xrange(m):
                slatec.dbspvn(t, k, k, 1, x[i], left[i], B[i], work, iwork)
        else:
            a = np.zeros((k, k))
            dbiatx = np.zeros((der+1, k))
            work = np.zeros((k+1)*(k+2)/2)
            for i in xrange(m):
                slatec.dbspvd(t, k, der+1, x[i], left[i], dbiatx.T, work)
                B[i] = dbiatx[der]
        return left, B

    def to_pp2(self):
        # this uses the PGS normalization
        t, c, k, p, n, d = self.spbrk()
        l = len(np.unique(t)) - 1

        b = np.zeros(l+1)
        scrtch = np.zeros((d, k, k))
        a = np.zeros((p, d, l, k))
        work = np.zeros(k*(n+3))

        c = np.ascontiguousarray(c.transpose((0, 2, 1)))
        for j in xrange(p):
            cj = c[j]
            aj = a[j]
            for i in xrange(d):
                slatec.dbsppp(t, cj[i], n, k, aj[i].T, b, l+1, work)
        a = np.ascontiguousarray(a.transpose((0, 2, 3, 1)))
        return PPSLA2(b, a)


class SplineND:
    def __init__(self, t, c):
        self.t, self.c = t, c

    def spval(self, x, y, derx=0, dery=0):
        t, c = self.t, self.c
        nx, ny = c.shape
        tx, ty = t
        kx, ky = tx.size - nx, ty.size - ny

        mx, my = x.size, y.size
        
        C = np.zeros((mx, ny))
        if derx == 0:
            pppack.spual(tx, c.reshape((1, nx, ny)).T, kx, x, C.reshape(1, mx, ny).T)
        else:
            pppack.spualder(tx, c.reshape((1, nx, ny)).T, kx, x, C.reshape(1, mx, ny).T, derx)

        Z = np.zeros((mx, my))
        if dery == 0:
            pppack.spual(ty, C.reshape(mx, ny, 1).T, ky, y, Z.reshape(mx, my, 1).T)
        else:
            pppack.spualder(ty, C.reshape(mx, ny, 1).T, ky, y, Z.reshape(mx, my, 1).T, dery)
        return Z


test1 = test2 = test3 = False
if __name__ == "__main__":
    test1 = True

if test1:
    p, n, d, k, m, der = 2, 10, 6, 4, 101, 2
    #p, n, d, k, m, der = 20, 1000, 20, 4, 10000, 1
    c = np.zeros((p, n, d))
    for i in xrange(d): c[:,n-1-i,i] = 1.
    #c = np.random.rand(p, n, d)

    #knots = np.arange(n-k+2.)
    knots = np.array((0., 1.2, 2.1, 2.1, 4.3, 5.2, 6.1, 7.))
    t = augknt(knots, k)
    sp = Spline.from_knots_coefs(t, c)

    x = np.linspace(knots[0], knots[-1], m)

    dsp = sp.deriv(der)

    y = dsp.spval(x)
    
    pp = sp.to_pp()
    dpp = pp.deriv(der)
    y2 = dpp.ppual(x)

    sp_pgs = SplineSLA.from_knots_coefs(t, c)

    y3 = sp_pgs.spval(x, der=der)
    y3b = sp_pgs.spval2(x, der=der)

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

    assert np.allclose(pp_pgs.a, pp_pgsb.a)

    #"""
    s = 0
    from matplotlib.pyplot import plot, show
    plot(x, np.c_[y[s], y2[s], y3[s], y4[s], y5[s], y6[s], y7[s], y8[s]], '.-')
    show()
    #"""

if test2:
    tx = np.array((0., 0., 0., 0., 2., 4., 4., 4., 4.))
    ty = np.array((0., 0., 0., 0., 2., 3., 5., 5., 5., 5.))

    """
    coefs = np.array(
       (( 1.00000,  1.06298, -0.25667, -1.40455, -0.42843,  0.28366),
        ( 1.07407,  1.14172, -0.27569, -1.50859, -0.46016,  0.30467),
        (-0.72984, -0.77581,  0.18733,  1.02510,  0.31269, -0.20703),
        (-1.27897, -1.35952,  0.32828,  1.79638,  0.54795, -0.36280),
        (-0.65364, -0.69481,  0.16777,  0.91808,  0.28004, -0.18541)))
    """
    coefs = np.random.rand(tx.size-kx, ty.size-ky)

    sp = SplineND((tx, ty), coefs)

    x = np.linspace(0., 4., 10)
    y = np.linspace(0., 5., 12)

    Z = sp.spval(x, y)

    from pytokamak.utils import splines
    sp2 = splines.Spline2D(data=dict(tck=(ty, tx, coefs.T.ravel()), degrees=(3, 3)))
    Z2 = sp2(y, x).T

    print (Z - Z2).ptp()

    from matplotlib.pyplot import figure, show
    from mpl_toolkits.mplot3d import Axes3D
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x[:,None], y[None,:], Z, rstride=1, cstride=1, cmap='jet')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    show()

if test3:
    nx, ny = 20, 25
    kx, ky = 3, 3
    tx = augknt(np.arange(nx-kx+2.), kx)
    ty = augknt(np.arange(ny-ky+2.), ky)

    c = np.random.rand(nx, ny)

    sp = SplineND((tx, ty), c)

    x = np.linspace(tx[0], tx[-1], 140)
    y = np.linspace(ty[0], ty[-1], 100)

    Z = sp.spval(x, y)

    from pytokamak.utils import splines
    sp2 = splines.Spline2D(data=dict(tck=(ty, tx, c.T.ravel()), degrees=(ky-1, kx-1)))
    Z2 = sp2(y, x).T

    print (Z - Z2).ptp()

    from scipy.ndimage import map_coordinates

    coordinates = np.meshgrid(x+0.5, y+0.5)

    Z3 = map_coordinates(c, coordinates, order=2, prefilter=False).T

    i = y.size / 2

    from matplotlib.pyplot import plot, show
    plot(x, Z[:,i], x, Z3[:,i])
    show()

