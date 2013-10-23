      subroutine spder( t, bcoef, n, k, d, jderiv )
      integer jderiv,n,k,d,   i,dd
      real t(n+k),bcoef(d,n),   h,fkmj

      do 50 j=1,jderiv
        fkmj = float(k-j)
        do 50 i=1,n-j
           h = fkmj/(t(i+k)-t(i+j))
           do 50 dd=1,d
              bcoef(dd,i) = h*(bcoef(dd,i+1)-bcoef(dd,i))
   50 continue
      end


