      subroutine spder( t, bcoef, n, k, d, l, jderiv )
      integer jderiv,n,k,d,l,   i,dd,ll
      real t(n+k),bcoef(d,n,l),   h,fkmj

      if (d .eq. 1 .and. l .eq. 1)         go to 99
      if (d .eq. 1)                        go to 199
      if (l .eq. 1)                        go to 299

      do 50 j=1,jderiv
        fkmj = float(k-j)
        do 50 i=1,n-j
          h = fkmj/(t(i+k)-t(i+j))
          do 50 ll=1,l
            do 50 dd=1,d
   50         bcoef(dd,i,ll) = h*(bcoef(dd,i+1,ll)-bcoef(dd,i,ll))
      return
   99 call spder0 (t, bcoef(1,:,1), n, k, jderiv)
      return
  199 call spderl (t, bcoef(1,:,:), n, k, l, jderiv)
      return
  299 call spderd (t, bcoef(:,:,1), n, k, d, jderiv)
      return
      end

      subroutine spder0( t, bcoef, n, k, jderiv )
      integer jderiv,n,k,   i
      real t(n+k),bcoef(n),   fkmj

      do 60 j=1,jderiv
        fkmj = float(k-j)
        do 60 i=1,n-j
   60     bcoef(i) = fkmj*(bcoef(i+1)-bcoef(i))/(t(i+k)-t(i+j))
      end

      subroutine spderl( t, bcoef, n, k, l, jderiv )
      integer jderiv,n,k,l,   i,ll
      real t(n+k),bcoef(n,l),   h,fkmj

      do 70 j=1,jderiv
        fkmj = float(k-j)
        do 70 i=1,n-j
          h = fkmj/(t(i+k)-t(i+j))
          do 70 ll=1,l
   70       bcoef(i,ll) = h*(bcoef(i+1,ll)-bcoef(i,ll))
      end

      subroutine spderd( t, bcoef, n, k, d, jderiv )
      integer jderiv,n,k,d,   i,dd
      real t(n+k),bcoef(d,n),   h,fkmj

      do 80 j=1,jderiv
        fkmj = float(k-j)
        do 80 i=1,n-j
          h = fkmj/(t(i+k)-t(i+j))
          do 80 dd=1,d
   80       bcoef(dd,i) = h*(bcoef(dd,i+1)-bcoef(dd,i))
      end

