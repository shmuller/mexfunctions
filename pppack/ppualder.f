      subroutine ppualder( break, coef, l, k, d, x, m, jderiv, y )
          integer jderiv,l,k,d,m,   i,dd,j,mm,ndummy,ifac
      real break(l+1),coef(d,k,l),x(m),y(d,m),   f,h
      
      if (jderiv .ge. k)            go to 99

      ifac = 1
      do 4 j=2,jderiv
    4   ifac = ifac*j

      do 50 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        do 8 dd=1,d
    8      y(dd,mm) = coef(dd,k,i)
        do 30 j=k-1,1+jderiv,-1
           f = h*j/(j-jderiv)
           do 10 dd=1,d
   10         y(dd,mm) = y(dd,mm)*f + coef(dd,j,i)
   30   continue
        if (ifac .eq. 1)            go to 50
        do 40 dd=1,d
   40      y(dd,mm) = y(dd,mm)*ifac 
   50 continue
      return
   99 do 199 mm=1,d*m
  199   y(mm,1) = 0.
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

