import numpy as np
diff, find, cat, zeros, ones = \
        np.diff, np.flatnonzero, np.concatenate, np.zeros, np.ones

def searchsorted(meshsites, sites):
    index = np.argsort(cat((meshsites, sites)))
    return find(index > len(meshsites)-1) - np.arange(1, len(sites)+1)

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

    def ppmak(self, brks, coefs, d):
        self.form = 'pp'
        self.brks = brks
        self.coefs = coefs
        self.pieces = len(brks) - 1
        self.order = coefs.shape[0]
        self.dim = d

    def ppbrk(self):
        return self.brks, self.coefs, self.pieces, self.order, self.dim

    def deriv(self, dorder=1):
        b, c, l, k, d = self.ppbrk()
        knew = k - dorder
        if knew > 0:
            fact = np.ones(knew)
            expo = np.arange(k-1, 0, -1)
            for j in xrange(dorder):
                fact[:] *= expo[j:j+knew]
            cnew = c[:knew] * fact[:,None]
        else:
            cnew = np.zeros((1, d*l))
        return self.__class__(b, cnew, d)

    def ppual(self, x):
        if any(diff(x) < 0):
            tosort = True
            ix = np.argsort(x)
            xs = x[ix]
        else:
            xs = x.copy()

        b, c, l, k, d = self.ppbrk()

        index = self.get_index(x)
        xs[:] -= b[index]

        xs = xs[:,None]
        c = c.reshape((k, l, d))
        v = c[0,index]
        for i in xrange(1, k):
            v = xs * v + c[i,index]
        return v

    def get_index(self, sites):
        ind = searchsorted(self.brks[:-1], sites)
        ind[ind < 0] = 0
        return ind

    def to_pp_pgs(self):
        b, c, l, k, d = self.ppbrk()
        coef = c[::-1].reshape((k,l,d)).transpose((1,0,2))
        return PPPGS(b, coef, d)


class Spline(object):
    def __init__(self, *args, **kw):
        self.spmak(*args, **kw)
    
    @classmethod
    def from_knots_coefs(cls, knots, coefs):
        self = cls.__new__(cls)
        self.spmak(knots, coefs)
        return self

    def spmak(self, knots, coefs, sizec=None):
        if sizec is None:
            sizec = coefs.shape
        m = 1
        if len(sizec) == m:
            sizec = sizec + (1,)

        sizeval = sizec[m:]
        sizec = sizec[:m] + (np.prod(sizeval),)
        coefs = coefs.reshape(sizec)

        knots, coefs, k, sizec = self.chckknt(knots, coefs, sizec)

        self.form = 'B-'
        self.knots = knots
        self.coefs = coefs
        self.number = sizec[0]
        self.order = k
        self.dim = sizeval

    def chckknt(self, knots, coefs, sizec):
        n = sizec[0]
        k = len(knots) - n
        
        # TODO some checks
        # TODO throw out trivial B-splines

        return knots, coefs, k, sizec

    def spbrk(self):
        return self.knots, self.coefs, self.number, self.order, self.dim
    
    def setup_tx_b(self, t, x, k, d, inter, backwd=False):
        if backwd:
            o = np.arange(k-1, 1-k, -1)
            kd = np.arange(0, -k*d, -d)
        else:
            o = np.arange(2-k, k)
            kd = np.arange(d*(1-k), d, d)

        tx = t[inter[None,:]+o[:,None]] - x[None,:]

        b = (d*(inter+1)-1)[None,:,None] + np.arange(1-d, 1)[None,None,:] \
          + kd[:,None,None]

        tx = tx[:,:,None].repeat(d, axis=2).reshape((2*(k-1),-1))
        b = b.reshape((k, -1))
        return tx, b

    def spval(self, x):
        t, a, n, k, dim = self.spbrk()
        d = prod(dim)

        index = searchsorted(t[:n], x)
        index[index < k-1] = k-1
        if k == 1:
            b = a[index]
        else:
            tx, b = self.setup_tx_b(t, x, k, d, index, backwd=False)
            b = a.ravel()[b]
            self.sprval(tx, b, k)
            b.resize((x.size,) + dim)
        return b

    def deriv(self, dorder=1):
        t, a, n, k, dim = self.spbrk()
        d = prod(dim)

        knew = k - dorder;
        if knew <= 0:
            t = t[k-1:n+1]
            a = zeros((n-k,) + dim)
        else:
            z = zeros((1, d))
            for j in xrange(k-1, knew-1, -1):
                tt = t[j+1:j+n] - t[1:n]
                i = find(tt > 0)
                temp = diff(a, axis=0)
                a = temp[i] * (j / tt[i])[:,None]
                t = cat((t[i+1], t[n+1:n+j+2]))
                n = len(i)
        return self.from_knots_coefs(t, a)

    def to_pp(self, backwd=True):
        t, a, n, k, dim = self.spbrk()
        d = prod(dim)

        inter = find(diff(t) > 0)
        if k == 1:
            b = a[inter].ravel()
        else:
            tx, b = self.setup_tx_b(t, t[inter], k, d, inter, backwd=backwd)
            b = a.ravel()[b]
            b = self.sprpp(tx, b, k, backwd=backwd)
        return PP(cat((t[inter], t[inter[-1]+1:inter[-1]+2])), b, d)

    def sprval(self, tx, b, k):
        for r in xrange(1, k):
            for i in xrange(k-r):
                tx0 = tx[i+r-1]
                tx1 = tx[i+k-1]
                b[i] = (tx1 * b[i] - tx0 * b[i+1]) / (tx1 - tx0)

    def sprval_bw(self, tx, b, k):
        for r in xrange(1, k):
            for j in xrange(k-1, r-1, -1):
                tx0 = tx[j-1+k-r]
                tx1 = tx[j-1]
                b[j] = (tx1 * b[j] - tx0 * b[j-1]) / (tx1 - tx0)

    def sprdif(self, tx, b, k):
        for r in xrange(1, k):
            factor = float(k-r) / r
            for i in xrange(k-1, r-1, -1):
                b[i] = (b[i] - b[i-1]) * factor / tx[i+k-1-r]

    def sprdif_bw(self, tx, b, k):
        for r in xrange(1, k):
            factor = float(k-r) / r
            for j in xrange(k-r):
                b[j] = (b[j] - b[j+1]) * factor / tx[j-1+r]

    def sprpp(self, tx, b, k, backwd=True):
        if backwd:
            self.sprval_bw(tx, b, k)
            self.sprdif_bw(tx, b, k)
            return b
        else:
            self.sprval(tx, b, k)
            self.sprdif(tx, b, k)
            return b[::-1]


import pppack

class PPPGS(PP):
    def ppmak(self, *args, **kw):
        PP.ppmak(self, *args, **kw)
        self.order = self.coefs.shape[1]

    def deriv(self, dorder=1):
        b, c, l, k, d = self.ppbrk()
        knew = k - dorder
        if knew > 0:
            fact = np.ones(knew)
            expo = np.arange(k-1, 0, -1)
            for j in xrange(dorder):
                fact[:] *= expo[j:j+knew]
            cnew = c[:,k-knew:] * fact[None,::-1,None]
        else:
            cnew = np.zeros((l, 1, d))
        return self.__class__(b, cnew, d)

    def ppual(self, x, der=0, fast=True):
        b, c, l, k, d = self.ppbrk()

        m = x.size
        y = np.zeros((m, d))
        if der == 0 and fast:
            pppack.ppual(b, c.T, x, y.T)
        else:
            pppack.ppualder(b, c.T, x, der, y.T)
        return y


class PPPGS2(PPPGS):
    def ppual(self, x, der=0):
        # this uses the PGS normalization
        b, c, l, k, d = self.ppbrk()

        m = x.size
        y = np.zeros((m, d))
        c = c.transpose((2,0,1)).copy()
        for i in xrange(m):
            yi = y[i]
            for j in xrange(d):
                yi[j] = pppack.ppvalu(b, c[j].T, x[i], der)
        return y


class SplinePGS(Spline):
    def spval(self, x, der=0):
        t, a, n, k, dim = self.spbrk()
        d = prod(dim)

        m = x.size
        y = np.zeros((m, d))
        pppack.spualder(t, a.reshape((1, n, d)).T, k, x, y.reshape(1, m, d).T, der)
        return y

    def deriv(self, dorder=1):
        t, a, n, k, dim = self.spbrk()
        d = prod(dim)
        if k <= dorder:
            t = t[k-1:n+1]
            anew = zeros((n-k,) + dim)
        else:
            anew = a.copy()
            pppack.spder(t, anew.reshape((1, n, d)).T, k, dorder)
            anew.resize((n-dorder,) + dim)
            t = t[dorder:n+k-dorder]
        return self.from_knots_coefs(t, anew)

    def to_pp(self):
        t, a, n, k, dim = self.spbrk()
        d = prod(dim)
        l = n+1-k

        b = np.zeros(l+1)
        scrtch = np.zeros((k, k, d))
        coef = np.zeros((l, k, d))

        l = pppack.bsplppd(t, a.reshape((1, n, d)).T, scrtch.reshape((k, k, 1, d)).T, 
                           b, coef.reshape((1, l, k, d)).T)
        if l < n+1-k:
            b.resize(l+1)
            coef.resize((l, k, d))
        return PPPGS(b, coef, d)

    def to_pp2(self):
        # this uses the PGS normalization
        t, a, n, k, dim = self.spbrk()
        d = prod(dim)
        l = n+1-k

        b = np.zeros(l+1)
        scrtch = np.zeros((d, k, k))
        coef = np.zeros((l, k, d))

        c = a.T.copy()
        pppack.bspp2d(t, c.T, scrtch.T, b, coef.T)
        l = len(np.unique(t)) - 1
        if l < n+1-k:
            b.resize(l+1)
            coef.resize((l, k, d))
        return PPPGS2(b, coef, d)


class SplineND:
    def __init__(self, knots, coefs):
        self.knots, self.coefs = knots, coefs

    def spval(self, x, y, derx=0, dery=0):
        t, c = self.knots, self.coefs
        nx, ny = c.shape
        tx, ty = t
        kx, ky = tx.size - nx, ty.size - ny

        mx, my = x.size, y.size
        
        C = np.zeros((mx, ny))
        pppack.spualder(tx, c.reshape((1, nx, ny)).T, kx, x, C.reshape(1, mx, ny).T, derx)

        Z = np.zeros((mx, my))
        pppack.spualder(ty, C.reshape(mx, ny, 1).T, ky, y, Z.reshape(mx, my, 1).T, dery)
        return Z


test1 = test2 = False
if __name__ == "__main__":
    test1 = True

if test1:
    n = 10
    #c = np.arange(2.*n).reshape(2, n).T.copy()
    m = 6
    c = np.zeros((n, m))
    for i in xrange(m): c[i,i] = 1.

    k = 4
    #knots = np.arange(n-k+2.)
    knots = np.array((0., 1.2, 2.5, 2.5, 4.3, 5.2, 6.1, 7.))
    t = augknt(knots, k)
    sp = Spline.from_knots_coefs(t, c)

    x = np.linspace(knots[0], knots[-1], 100)

    der = 1
    dsp = sp.deriv(der)

    y = dsp.spval(x)
    
    pp = sp.to_pp()
    y2 = pp.deriv(der).ppual(x)

    sp_pgs = SplinePGS.from_knots_coefs(t, c)

    y3 = sp_pgs.spval(x, der=der)

    dsp_pgs = sp_pgs.deriv(der)

    y4 = dsp_pgs.spval(x)

    pp_pgs = sp_pgs.to_pp()
    y5 = pp_pgs.ppual(x, der=der, fast=False)

    dpp_pgs = pp_pgs.deriv(der)
    y6 = dpp_pgs.ppual(x, fast=False)

    pp_pgs2 = sp_pgs.to_pp2()
    y7 = pp_pgs2.ppual(x, der=der)

    pp_pgsb = pp.to_pp_pgs()
    y8 = pp_pgsb.ppual(x, der=der)

    assert np.allclose(pp_pgs.coefs, pp_pgsb.coefs)

    from matplotlib.pyplot import plot, show
    plot(x, np.c_[y, y2, y3, y4, y5, y6, y7, y8], '.-')
    show()

if test2:
    tx = np.array((0., 0., 0., 0., 2., 4., 4., 4., 4.))
    ty = np.array((0., 0., 0., 0., 2., 3., 5., 5., 5., 5.))

    coefs = np.array(
       (( 1.00000,  1.06298, -0.25667, -1.40455, -0.42843,  0.28366),
        ( 1.07407,  1.14172, -0.27569, -1.50859, -0.46016,  0.30467),
        (-0.72984, -0.77581,  0.18733,  1.02510,  0.31269, -0.20703),
        (-1.27897, -1.35952,  0.32828,  1.79638,  0.54795, -0.36280),
        (-0.65364, -0.69481,  0.16777,  0.91808,  0.28004, -0.18541)))

    sp = SplineND((tx, ty), coefs)

    x = np.linspace(0., 4., 10)
    y = np.linspace(0., 5., 12)

    Z = sp.spval(x, y)

    from matplotlib.pyplot import figure, show
    from mpl_toolkits.mplot3d import Axes3D
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x[:,None], y[None,:], Z, rstride=1, cstride=1, cmap='jet')
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    show()





