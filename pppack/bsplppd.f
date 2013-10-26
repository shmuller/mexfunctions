      subroutine bsplppd (t, bcoef, n, k, d, p, scrtch, break, coef, l)
cf2py intent(out) l
c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
calls  bsplvb
c
converts the b-representation  t, bcoef, n, k  of some spline into its
c  pp-representation  break, coef, l, k .
c
c******  i n p u t  ******
c  t.....knot sequence, of length  n+k
c  bcoef.....b-spline coefficient sequence, of length  n
c  n.....length of  bcoef  and  dimension of spline space  spline(k,t)
c  k.....order of the spline
c
c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
c        sed by the arbitrary dimension statement for  biatx  below, but
c        is  n o w h e r e   c h e c k e d   for.
c
c******  w o r k   a r e a  ******
c  scrtch......of size  (k,k) , needed to contain bcoeffs of a piece of
c        the spline and its  k-1  derivatives
c
c******  o u t p u t  ******
c  break.....breakpoint sequence, of length  l+1, contains (in increas-
c        ing order) the distinct points in the sequence  t(k),...,t(n+1)
c  coef.....array of size (k,l), with  coef(i,j) = (i-1)st derivative of
c        spline at break(j) from the right
c  l.....number of polynomial pieces which make up the spline in the in-
c        terval  (t(k), t(n+1))
c
c******  m e t h o d  ******
c     for each breakpoint interval, the  k  relevant b-coeffs of the
c  spline are found and then differenced repeatedly to get the b-coeffs
c  of all the derivatives of the spline on that interval. the spline and
c  its first  k-1  derivatives are then evaluated at the left end point
c  of that interval, using  bsplvb  repeatedly to obtain the values of
c  all b-splines of the appropriate order at that point.
c
      integer k,l,n,d,p,   i,jp1,kmax,kmj,left,dd,pp,dp
      parameter (kmax = 20)
      real t(n+k),bcoef(d,n,p),break(n+2-k),coef(d,k,n+1-k,p),
     *     scrtch(d,p,k,k),biatx(kmax),s
c
      dp = d*p
      l = 0
      break(1) = t(k)
      do 50 left=k,n
c                                find the next nontrivial knot interval.
         if (t(left+1) .eq. t(left))    go to 50
         l = l + 1
         break(l+1) = t(left+1)
         if (k .gt. 1)                  go to 9
         do 5 pp=1,p
           do 5 dd=1,d
    5        coef(dd,1,l,pp) = bcoef(dd,left,pp)
                                        go to 50
c        store the k b-spline coeff.s relevant to current knot interval
c                             in  scrtch(.,1) .
    9    do 10 i=1,k
           do 10 pp=1,p
             do 10 dd=1,d
   10          scrtch(dd,pp,i,1) = bcoef(dd,left-k+i,pp)
c
c        for j=1,...,k-1, compute the  k-j  b-spline coeff.s relevant to
c        current knot interval for the j-th derivative by differencing
c        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
         call calcbspl ( t(left+1), scrtch, k, dp )
c
c        for  j = 0, ..., k-1, find the values at  t(left)  of the  j+1
c        b-splines of order  j+1  whose support contains the current
c        knot interval from those of order  j  (in  biatx ), then comb-
c        ine with the b-spline coeff.s (in scrtch(.,k-j) ) found earlier
c        to compute the (k-j-1)st derivative at  t(left)  of the given
c        spline.
c           note. if the repeated calls to  bsplvb  are thought to gene-
c        rate too much overhead, then replace the first call by
c           biatx(1) = 1.
c        and the subsequent call by the statement
c           j = jp1 - 1
c        followed by a direct copy of the lines
c           deltar(j) = t(left+j) - x
c                  ......
c           biatx(j+1) = saved
c        from  bsplvb . deltal(kmax)  and  deltar(kmax)  would have to
c        appear in a dimension statement, of course.
c
         call bsplvb ( t, 1, 1, t(left), left, biatx )
         do 25 pp=1,p
           do 25 dd=1,d
   25        coef(dd,k,l,pp) = scrtch(dd,pp,1,k)
         do 30 jp1=2,k
            call bsplvb ( t, jp1, 2, t(left), left, biatx )
            kmj = k+1 - jp1
            do 30 pp=1,p
              do 30 dd=1,d
                s = 0.
                do 28 i=1,jp1
   28              s = s + biatx(i)*scrtch(dd,pp,i,kmj)
   30           coef(dd,kmj,l,pp) = s
   50 continue
      end

      subroutine calcbspl ( t, scr, k, d )
      integer k,d,   jp1,j,i,kmj,dd
      real t(1),scr(d,k,k),   diff,fact,h

      do 100 jp1=2,k
        j = jp1-1
        kmj = k-j
        fact = float(kmj)/j
        do 100 i=1,kmj
          diff = t(i)-t(i-kmj)
          if (diff .le. 0.)             go to 100
          h = fact/diff
          do 90 dd=1,d
   90       scr(dd,i,jp1) = (scr(dd,i+1,j)-scr(dd,i,j))*h
  100 continue
      end
