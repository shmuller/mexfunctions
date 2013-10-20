      subroutine ppualder( break, coef, l, k, d, x, m, jderiv, y )
      integer jderiv,l,k,d,m,   i,dd,j,mm,ndummy,n,kmax,ifac(10)
      real break(l+1),coef(d,k,l),x(m),y(d,m),   f,h,hf
      parameter(kmax=10)
      logical first
      data first/.true./
      save first,ifac
      if (first) then
          ifac(1) = 1
          ifac(2) = 1
          do 3 j=3,kmax
    3        ifac(j) = ifac(j-1)*(j-1)
          first = .false.
      end if

      do 6 mm=1,d*m
    6    y(mm,1) = 0.

      if (jderiv .ge. k) return

      do 50 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        f = k - jderiv
        do 50 j=k,1+jderiv,-1
           hf = h/f
           do 10 dd=1,d
   10         y(dd,mm) = y(dd,mm)*hf + coef(dd,j,i)*ifac(j)
           f = f - 1.
   50 continue
      end
