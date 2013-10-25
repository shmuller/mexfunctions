      subroutine spualder( t, bcoef, n, k, d, l, x, m, y, jderiv )
      integer n,k,d,l,m,   dd,ll,mm
      real t(n+k),bcoef(d,n,l),x(m),y(d,m,l),  xm

      do 50 mm=1,m
        xm = x(mm)
        do 50 ll=1,l
          do 50 dd=1,d
             y(dd,mm,ll) = bvalue (t, bcoef(dd,:,ll), n, k, xm, jderiv)
   50 continue
      end
