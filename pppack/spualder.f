      subroutine spualder( t, bcoef, n, k, d, x, m, y, jderiv )
      integer n,k,d,m,   dd,mm
      real t(n+k),bcoef(d,n),x(m),y(d,m)

      do 50 mm=1,m
        do 50 dd=1,d
           y(dd,mm) = bvalue ( t, bcoef(dd,:), n, k, x(mm), jderiv )
   50 continue
      end
