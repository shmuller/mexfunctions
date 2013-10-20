      subroutine ppualder( break, coef, l, k, d, x, m, jderiv, y )
      integer jderiv,l,k,d,m,   i,dd,j,mm,ndummy,n
      real break(l+1),coef(d,k,l),x(m),y(d,m),   f,h,factor

      do 1 mm=1,d*m
    1    y(mm,1) = 0.

      if (jderiv .ge. k) return

      do 50 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        f = k - jderiv
        do 30 j=k,1+jderiv,-1
           factor = 1.
           do 7 n=1,j-1
    7         factor = factor*n
           do 10 dd=1,d
   10         y(dd,mm) = (y(dd,mm)/f)*h + coef(dd,j,i)*factor
           f = f - 1.
   30   continue
   50 continue
      end
