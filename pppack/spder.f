      subroutine spder( t, bcoef, n, k, d, l, jderiv )
      integer jderiv,n,k,d,l,   i,dd,ll
      real t(n+k),bcoef(d,n,l),   h,fkmj

      do 50 j=1,jderiv
        fkmj = float(k-j)
        do 50 i=1,n-j
          h = fkmj/(t(i+k)-t(i+j))
          do 50 ll=1,l
            do 50 dd=1,d
   50         bcoef(dd,i,ll) = h*(bcoef(dd,i+1,ll)-bcoef(dd,i,ll))
      end


