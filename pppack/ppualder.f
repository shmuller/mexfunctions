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


c   y = 0
c   y = y/(k-0-j)*x + a*(k-1)*(k-2)*(k-3)
c   y = y/(k-1-j)*x + b*(k-2)*(k-3) 
c   y = y/(k-2-j)*x + c*(k-3)
c   y = y/(k-3-j)*x + d

c   y = ((a*(k-1)*(k-2)*(k-3)/(k-1-j)*x + b*(k-2)*(k-3))/(k-2-j)*x
c     + c*(k-3))/(k-3-j)*x + d

c   y = ((a*(k-1)*(k-2)*(k-3)/(k-1-j)/(k-2-j)/(k-3-j)*x +
c         b*(k-2)*(k-3)/(k-2-j)/(k-3-j))*x +
c         c*(k-3)/(k-3-j))*x +
c         d

c   y = ((a*r1*r2*r3*x + b*r2*r3)*x + c*r3)*x + d
c
c   ri := (k-i)/(k-i-j)

c   y = ((a*r1*x + b)*r2*x +c)*r3*x + d
c
c
c

