!DECK DBSPED
      SUBROUTINE DBSPED (T, AD, N, K, NDERIV, X, INEV, SVALUE, WORK)
!***BEGIN PROLOGUE  DBSPED
!***PURPOSE  Calculate the value of the spline and its derivatives from
!            the B-representation.
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      DOUBLE PRECISION (BSPEV-S, DBSPED-D)
!***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract    **** a double precision routine ****
!         DBSPED is the BSPLEV routine of the reference.
!
!         DBSPED calculates the value of the spline and its derivatives
!         at X from the B-representation (T,A,N,K) and returns them in
!         SVALUE(I),I=1,NDERIV, T(K) .LE. X .LE. T(N+1).  AD(I) can be
!         the B-spline coefficients A(I), I=1,N) if NDERIV=1.  Otherwise
!         AD must be computed before hand by a call to DBSPDR (T,A,N,K,
!         NDERIV,AD).  If X=T(I),I=K,N), right limiting values are
!         obtained.
!
!         To compute left derivatives or left limiting values at a
!         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!         DBSPED calls DINTRV, DBSPVN
!
!     Description of Arguments
!
!         Input      T,AD,X, are double precision
!          T       - knot vector of length N+K
!          AD      - vector of length (2*N-NDERIV+1)*NDERIV/2 containing
!                    the difference table from DBSPDR.
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K .GE. 1
!          NDERIV  - number of derivatives, 1 .LE. NDERIV .LE. K.
!                    NDERIV=1 gives the zero-th derivative =
!                    function value
!          X       - argument, T(K) .LE. X .LE. T(N+1)
!          INEV    - an initialization parameter which must be set
!                    to 1 the first time DBSPED is called.
!
!         Output     SVALUE,WORK are double precision
!          INEV    - INEV contains information for efficient process-
!                    ing after the initial call and INEV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INEV parameters.
!          SVALUE  - vector of length NDERIV containing the spline
!                    value in SVALUE(1) and the NDERIV-1 derivatives
!                    in the remaining components.
!          WORK    - work vector of length 3*K
!
!     Error Conditions
!         Improper input is a fatal error.
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  DBSPVN, DINTRV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSPED
!
      INTEGER I,ID,INEV,IWORK,JJ,K,KP1,KP1MN,L,LEFT,LL,MFLAG, &
       N, NDERIV
      DOUBLE PRECISION AD, SVALUE, SUM, T, WORK, X
!     DIMENSION T(N+K)
      DIMENSION T(*), AD(*), SVALUE(*), WORK(*), INEV(*)
!***FIRST EXECUTABLE STATEMENT  DBSPED
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      IF(NDERIV.LT.1 .OR. NDERIV.GT.K) GO TO 115
      ID = NDERIV
      CALL DINTRV(T, N+1, X, INEV(1), I, MFLAG)
      IF (X.LT.T(K)) GO TO 110
      IF (MFLAG.EQ.0) GO TO 30
      IF (X.GT.T(I)) GO TO 110
   20 IF (I.EQ.K) GO TO 120
      I = I - 1
      IF (X.EQ.T(I)) GO TO 20
!
! *I* HAS BEEN FOUND IN (K,N) SO THAT T(I) .LE. X .LT. T(I+1)
!     (OR .LE. T(I+1), IF T(I) .LT. T(I+1) = T(N+1) ).
   30 KP1 = K + 1
      KP1MN = KP1 - ID
      CALL DBSPVN(T, KP1MN, K, 1, X, I, WORK(1),WORK(KP1),IWORK)
!     ADIF(LEFTPL,ID) = AD(LEFTPL-ID+1 + (2*N-ID+2)*(ID-1)/2)
!     LEFTPL = LEFT + L
      SUM = 0.0D0
      LL = I - K + 1 + (N+N-ID+2)*(ID-1)/2
      DO 50 L=1,KP1MN
        SUM = SUM + WORK(L)*AD(LL)
        LL = LL + 1
   50 CONTINUE
      SVALUE(1) = SUM
   60 RETURN


  100 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPED', 'K DOES NOT SATISFY K.GE.1', 2, &
         1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPED', 'N DOES NOT SATISFY N.GE.K', 2, &
         1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPED', &
         'X IS NOT IN T(K).LE.X.LE.T(N+1)', 2, 1)
      RETURN
  115 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPED', &
         'NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K', 2, 1)
      RETURN
  120 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPED', &
         'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
      RETURN
      END
