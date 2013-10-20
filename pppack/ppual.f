      subroutine ppual( break, coef, l, k, d, x, m, y )
      integer l,k,d,m,   i,dd,kk,mm,ndummy
      real break(l+1),coef(d,k,l),x(m),y(d,m),   h

      do 50 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        do 5 dd=1,d
    5      y(dd,mm) = coef(dd,k,i)
        do 10 kk=k-1,1,-1
           do 10 dd=1,d
   10         y(dd,mm) = y(dd,mm)*h + coef(dd,kk,i)
   50 continue
      end
