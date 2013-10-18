import numpy as np
diff, find, cat = np.diff, np.flatnonzero, np.concatenate

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
    
    def to_pp(self, backwd=False):
        t, a, n, k, d = self.spbrk()
        d = prod(d)
        m = 1

        index = find(diff(t) > 0); addl = k-index[0]-1; addr = index[-1]+1-n
        if addl > 0 or addr > 0:
            tt = np.zeros(addl + t.size + addr)
            tt[:addl] = t[0]
            tt[addl:-addr] = t
            tt[-addr:] = t[-1]
            t = tt

        inter = find(diff(t) > 0); l = len(inter)
        if k == 1:
            b = a[inter].ravel()
        else:
            if backwd:
                o = np.arange(k-1, 1-k, -1)
                kd = np.arange(0, -k*d, -d)
            else:
                o = np.arange(2-k, k)
                kd = np.arange(d*(1-k), d, d)

            i = inter[None,:]
            tx = t[i+o[:,None]] - t[i]

            b = (d*(inter+1)-1)[None,:,None] + np.arange(1-d, 1)[None,None,:] \
              + kd[:,None,None]

            tx = tx[:,:,None].repeat(d, axis=2).reshape((2*(k-1),-1))
            b = b.reshape((k, -1))

            b = a.ravel()[b]
            if backwd:
                b = self.sprpp_bw(tx, b)
            else:
                b = self.sprpp(tx, b)

        return PP(cat((t[inter], t[inter[-1]+1:inter[-1]+2])), b, d)

    def sprpp(self, tx, b):
        k = b.shape[0]
        for r in xrange(1, k):
            for i in xrange(k-r):
                tx0 = tx[i+r-1]
                tx1 = tx[i+k-1]
                b[i] = (tx1 * b[i] - tx0 * b[i+1]) / (tx1 - tx0)

        for r in xrange(1, k):
            factor = float(k-r) / r
            for i in xrange(k-1, r-1, -1):
                b[i] = (b[i] - b[i-1])*factor / tx[i+k-1-r]
        return b[::-1]

    def sprpp_bw(self, tx, b):
        k = b.shape[0]
        for r in xrange(1, k):
            for j in xrange(k-1, r-1, -1):
                tx0 = tx[j-1+k-r]
                tx1 = tx[j-1]
                b[j] = (tx1 * b[j] - tx0 * b[j-1]) / (tx1 - tx0)
        
        for r in xrange(1, k):
            factor = float(k-r) / r
            for j in xrange(k-r):
                b[j] = (b[j] - b[j+1]) * factor / tx[j-1+r]
        return b


if __name__ == "__main__":
    #c = np.arange(2.*n).reshape(2, n).T.copy()
    n = 10
    c = np.zeros((n, n))
    for i in xrange(n): c[i,i] = 1.

    k = 4
    knots = np.arange(n-k+2)
    sp = Spline.from_knots_coefs(augknt(knots, k), c)

    pp = sp.to_pp()

    x = np.linspace(knots[0], knots[-1], 100)

    y = pp.ppual(x)

