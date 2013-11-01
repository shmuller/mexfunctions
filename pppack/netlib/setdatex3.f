      subroutine setdatex3 ( icount ,      
     *                       ntau, tau, gtau, weight, totalw,
     *                       break, coef, l, k )
cf2py intent(out) ntau, totalw, l, k
c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
c  to be called in main program  l 2 m a i n .
c     this routine is set up to provide the specific data for example 3
c     in chapter xiv.
c
      integer icount,  i,k,l,ntau,lpkmax,ntmax,ltkmax
      real break,coef,gtau,step,tau,totalw,weight,round,x
      parameter (lpkmax=100,ntmax=200,ltkmax=2000)
      dimension tau(ntmax),gtau(ntmax),weight(ntmax)
      dimension break(lpkmax),coef(ltkmax)
C     common / data / ntau, tau(200),gtau(200),weight(200),totalw
C     common /approx/ break(100),coef(2000),l,k
      round(x) = float(ifix(x*100.))/100.
      if (icount .gt. 0)                stop
      icount = icount + 1
      ntau = 65
      step = 3./float(ntau-1)
      do 10 i=1,ntau
         tau(i) = (i-1)*step
         gtau(i) = round(exp(tau(i)))
   10    weight(i) = 1.
      totalw = ntau
      l = 1
      break(1) = tau(1)
      break(2) = tau(ntau)
      k = 3
                                        return
      end
