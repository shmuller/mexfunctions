      subroutine ppualder( break, coef, l, k, d, x, m, jderiv, y )
      integer jderiv,l,k,d,m,   i,dd,j,mm,ndummy
      real break(l+1),coef(d,k,l),x(m),y(d,m),   f,h

      if (jderiv .ge. k)               go to 99

      do 50 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        do 5 dd=1,d
    5      y(dd,mm) = coef(dd,k,i)
        do 50 j=k-1,1+jderiv,-1
           f = h*j/(j-jderiv)
           do 10 dd=1,d
   10         y(dd,mm) = y(dd,mm)*f + coef(dd,j,i)
   50 continue
      return
   99 do 199 mm=1,d*m
  199   y(mm,1) = 0.
      return
      end
