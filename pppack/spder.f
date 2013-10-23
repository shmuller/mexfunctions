      subroutine spder( t, bcoef, n, k, d, jderiv, dcoef)
          integer jderiv,n,k,d,   i,j
      real t(n+k),bcoef(d,n),dcoef(d,n-jderiv),   h

      do 50 i=1,n-jderiv
        h = (k-1)/(t(i+k)-t(i+1))
        do 50 j=1,d
            dcoef(j,i) = h*(bcoef(j,i+1)-bcoef(j,i))
   50 continue
      end


