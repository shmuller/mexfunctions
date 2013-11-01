      subroutine spual( t, bcoef, n, k, d, p, x, m, y )
      integer n,k,d,p,m,jderiv,   dd,pp,mm
      real t(n+k),bcoef(d,n,p),x(m),y(d,m,p),bvalue

      if (d .eq. 1 .and. p .eq. 1)         go to 99

      do 50 mm=1,m
   50   call bual (t, bcoef, n, k, d, p, x(mm), y(:,mm,:))
      return
   99 call eval0 (t, bcoef(1,:,1), n, k, x, m, y(1,:,1))
      return
      end

      subroutine eval0( t, bcoef, n, k, x, m, y )
      integer n,k,m,   mm
      real t(n+k),bcoef(n),x(m),y(m),bvalue
      do 60 mm=1,m
   60   y(mm) = bvalue (t, bcoef, n, k, x(mm), 0)
      end

