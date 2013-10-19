      subroutine bsplppd ( t, bcoef, n, k, d, scrtch, break, coef, l )
cf2py intent(out) l
      integer n,k,d,l,   i
      real t(n+k),bcoef(d,n),scrtch(k,k),break(n+2-k),coef(d,k,n+1-k)
      do 10 i=1,d
   10   call bsplpp(t, bcoef(i,:), n, k, scrtch, break, coef(i,:,:), l)
      return
      end
