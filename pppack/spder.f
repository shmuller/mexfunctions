      subroutine spder( t, bcoef, n, k, d, p, jderiv )
      integer jderiv,n,k,d,p,   i,j,dd,pp
      real t(n+k),bcoef(d,n,p),   h,fkmj

      if (d .eq. 1 .and. p .eq. 1)         go to 99
      if (d .eq. 1)                        go to 199
      if (p .eq. 1)                        go to 299

      do 50 j=1,jderiv
        fkmj = float(k-j)
        do 50 i=1,n-j
          h = fkmj/(t(i+k)-t(i+j))
          do 50 pp=1,p
            do 50 dd=1,d
   50         bcoef(dd,i,pp) = h*(bcoef(dd,i+1,pp)-bcoef(dd,i,pp))
      return
   99 call spder0 (t, bcoef(1,:,1), n, k, jderiv)
      return
  199 call spderp (t, bcoef(1,:,:), n, k, p, jderiv)
      return
  299 call spderd (t, bcoef(:,:,1), n, k, d, jderiv)
      return
      end

      subroutine spder0( t, bcoef, n, k, jderiv )
      integer jderiv,n,k,   i,j
      real t(n+k),bcoef(n),   fkmj

      do 60 j=1,jderiv
        fkmj = float(k-j)
        do 60 i=1,n-j
   60     bcoef(i) = fkmj*(bcoef(i+1)-bcoef(i))/(t(i+k)-t(i+j))
      end

      subroutine spderp( t, bcoef, n, k, p, jderiv )
      integer jderiv,n,k,p,   i,j,pp
      real t(n+k),bcoef(n,p),   h,fkmj

      do 70 j=1,jderiv
        fkmj = float(k-j)
        do 70 i=1,n-j
          h = fkmj/(t(i+k)-t(i+j))
          do 70 pp=1,p
   70       bcoef(i,pp) = h*(bcoef(i+1,pp)-bcoef(i,pp))
      end

      subroutine spderd( t, bcoef, n, k, d, jderiv )
      integer jderiv,n,k,d,   i,j,dd
      real t(n+k),bcoef(d,n),   h,fkmj

      do 80 j=1,jderiv
        fkmj = float(k-j)
        do 80 i=1,n-j
          h = fkmj/(t(i+k)-t(i+j))
          do 80 dd=1,d
   80       bcoef(dd,i) = h*(bcoef(dd,i+1)-bcoef(dd,i))
      end

