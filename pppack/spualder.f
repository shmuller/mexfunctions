      subroutine spualder( t, bcoef, n, k, d, l, x, m, y, jderiv )
      integer n,k,d,l,m,jderiv,   dd,ll,mm
      real t(n+k),bcoef(d,n,l),x(m),y(d,m,l),  xm

      if (d .eq. 1 .and. l .eq. 1)         go to 99
      if (d .eq. 1)                        go to 199
      if (l .eq. 1)                        go to 299

      do 50 mm=1,m
        xm = x(mm)
        do 50 ll=1,l
          do 50 dd=1,d
   50        y(dd,mm,ll) = bvalue (t, bcoef(dd,:,ll), n, k, xm, jderiv)
      return
   99 call eval0 (t, bcoef(1,:,1), n, k, x, m, y(1,:,1), jderiv)
      return
  199 call evall (t, bcoef(1,:,:), n, k, l, x, m, y(1,:,:), jderiv)
      return
  299 call evald (t, bcoef(:,:,1), n, k, d, x, m, y(:,:,1), jderiv)
      return
      end

      subroutine eval0( t, bcoef, n, k, x, m, y, jderiv )
      integer n,k,m,jderiv,   mm
      real t(n+k),bcoef(n),x(m),y(m)
      do 60 mm=1,m
   60   y(mm) = bvalue (t, bcoef, n, k, x(mm), jderiv)
      end

      subroutine evall( t, bcoef, n, k, l, x, m, y, jderiv )
      integer n,k,d,m,jderiv,   dd,mm
      real t(n+k),bcoef(n,l),x(m),y(m,l),   xm
      do 70 mm=1,m
        xm = x(mm)
        do 70 ll=1,l
   70     y(mm,ll) = bvalue (t, bcoef(:,ll), n, k, xm, jderiv)
      end

      subroutine evald( t, bcoef, n, k, d, x, m, y, jderiv )
      integer n,k,d,m,jderiv,   dd,mm
      real t(n+k),bcoef(d,n),x(m),y(d,m),   xm
      do 80 mm=1,m
        xm = x(mm)
        do 80 dd=1,d
   80     y(dd,mm) = bvalue (t, bcoef(dd,:), n, k, xm, jderiv)
      end

