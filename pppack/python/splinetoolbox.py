import numpy as np

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


class Spline:
    def __init__(self, *args, **kw):
        self.spmak(*args, **kw)
        
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
    
    def to_pp(self):
        diff, find = np.diff, np.flatnonzero
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
            c = a[inter].ravel()
        else:
            i = inter[None,:]
            o = np.arange(k-1, 1-k, -1)[:,None]
            tx = t[i+o] - t[i]

            b = (d*(inter+1)-1)[None,:,None] + np.arange(1-d, 1)[None,None,:] \
              + (d*np.arange(0, -k, -1))[:,None,None]

            tx = tx[:,:,None].repeat(d, axis=2).reshape((2*(k-1),-1))
            b = b.reshape((l, -1))

            b = a.ravel()[b]
            c = self.sprpp(tx, b)

    def sprpp(self, tx, b):
        k = b.shape[0]
        for r in xrange(1, k):
            for i in xrange(k-r):
                tx0 = tx[2*(k-1)-1-(i+r-1)]
                tx1 = tx[2*(k-1)-1-(i+k-1)]
                b[k-1-i] = (tx1 * b[k-1-i] - tx0 * b[k-1-(i+1)]) / (tx1 - tx0)
        
        for r in xrange(2, k+1):
            factor = float(k-(r-1))/float(r-1)
            for i in xrange(k-1, r-2, -1):
                b[k-1-i] = (b[k-1-i] - b[k-1-(i-1)])*factor / tx[2*(k-1)-1-(i+k-r)]

        print b.T
    


if __name__ == "__main__":
    c = np.arange(14.).reshape(2, 7).T.copy()

    sp = Spline(augknt(np.arange(5.), 4), c)

    pp = sp.to_pp()


