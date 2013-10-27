      subroutine ppual ( break, coef, l, k, d, p, x, m, y )
      integer l,k,d,p,m,   i,dd,pp,kk,mm,ndummy
      real break(l+1),coef(d,k,l,p),x(m),y(d,m,p),   h,s

      if (d .eq. 1 .and. p .eq. 1)         go to 99
      if (d .eq. 1)                        go to 199
      if (p .eq. 1)                        go to 299

      do 20 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        do 20 pp=1,p
          do 20 dd=1,d
            s = coef(dd,k,i,pp)
            do 10 kk=k-1,1,-1
   10         s = s*h + coef(dd,kk,i,pp)
   20       y(dd,mm,pp) = s
      return
   99 call ppual0 (break, coef(1,:,:,1), l, k, x, m, y(1,:,1))
      return
  199 call ppualp (break, coef(1,:,:,:), n, k, p, x, m, y(1,:,:))
      return
  299 call ppuald (break, coef(:,:,:,1), n, k, d, x, m, y(:,:,1))
      return
      end

      subroutine ppual0 ( break, coef, l, k, x, m, y )
      integer l,k,m,   i,kk,mm,ndummy
      real break(l+1),coef(k,l),x(m),y(m),   h,s

      do 40 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        s = coef(k,i)
        do 30 kk=k-1,1,-1
   30     s = s*h + coef(kk,i)
   40   y(mm) = s
      end

      subroutine ppualp ( break, coef, l, k, p, x, m, y )
      integer l,k,p,m,   i,pp,kk,mm,ndummy
      real break(l+1),coef(k,l,p),x(m),y(m,p),   h,s

      do 60 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        do 60 pp=1,p
          s = coef(k,i,pp)
          do 50 kk=k-1,1,-1
   50       s = s*h + coef(kk,i,pp)
   60     y(mm,pp) = s
      end

      subroutine ppuald ( break, coef, l, k, d, x, m, y )
      integer l,k,d,m,   i,dd,kk,mm,ndummy
      real break(l+1),coef(d,k,l),x(m),y(d,m),   h,s

      do 80 mm=1,m
        call interv ( break, l+1, x(mm), i, ndummy )
        h = x(mm) - break(i)
        do 80 dd=1,d
          s = coef(dd,k,i)
          do 70 kk=k-1,1,-1
   70       s = s*h + coef(dd,kk,i)
   80     y(dd,mm) = s
      end

