      subroutine ppual( break, coef, l, k, d, x, m, jderiv, y )
      integer jderiv,l,k,d,m,   i,dd,kk,mm,ndummy
      real break(l+1),coef(d,k,l),x(m),y(d,m),   fmmjdr,f,h
      fmmjdr = k - jderiv

      if (fmmjdr .le. 0.)               go to 99

      do 50 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        f = fmmjdr
        kk = k
    9   do 10 dd=1,d
   10      y(dd,mm) = (y(dd,mm)/f)*h + coef(dd,kk,i)
        kk = kk - 1
        f = f - 1.
   50 if (f .gt. 0.)                    go to 9
      return
   99 do 199 mm=1,d*m
  199   y(mm,1) = 0.
      return
      end
