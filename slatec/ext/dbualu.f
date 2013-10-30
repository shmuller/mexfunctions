*DECK DBUALU
      DOUBLE PRECISION FUNCTION DBUALU (T, A, N, K, IDERIV, X, INBV,
     +   WORK)
C***BEGIN PROLOGUE  DBUALU
C***PURPOSE  Evaluate the B-representation of a B-spline at X for the
C            function value or any of its derivatives.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (BVALU-S, DBUALU-D)
C***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract   **** a double precision routine ****
C         DBUALU is the BVALUE function of the reference.
C
C         DBUALU evaluates the B-representation (T,A,N,K) of a B-spline
C         at X for the function value on IDERIV=0 or any of its
C         derivatives on IDERIV=1,2,...,K-1.  Right limiting values
C         (right derivatives) are returned except at the right end
C         point X=T(N+1) where left limiting values are computed.  The
C         spline is defined on T(K) .LE. X .LE. T(N+1).  DBUALU returns
C         a fatal error message when X is outside of this interval.
C
C         To compute left derivatives or left limiting values at a
C         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
C
C         DBUALU calls DINTRV
C
C     Description of Arguments
C
C         Input      T,A,X are double precision
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
C                    IDERIV = 0 returns the B-spline value
C          X       - argument, T(K) .LE. X .LE. T(N+1)
C          INBV    - an initialization parameter which must be set
C                    to 1 the first time DBUALU is called.
C
C         Output     WORK,DBUALU are double precision
C          INBV    - INBV contains information for efficient process-
C                    ing after the initial call and INBV must not
C                    be changed by the user.  Distinct splines require
C                    distinct INBV parameters.
C          WORK    - work vector of length 3*K.
C          DBUALU  - value of the IDERIV-th derivative at X
C
C     Error Conditions
C         An improper input is a fatal error
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  DINTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBUALU
C
      INTEGER I, IDERIV, INBV, IP1, K, IP1MK, KMIDER, KM1, MFLAG, N
      DOUBLE PRECISION A, FKMJ, T, WORK, X
      DIMENSION T(*), A(*), WORK(*)
C***FIRST EXECUTABLE STATEMENT  DBUALU
      DBUALU = 0.0D0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
C
C *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
C
C *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
C     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
C
   20 IP1 = I + 1
      IP1MK = IP1 - K
      CALL DINIT (A(IP1MK), K, WORK)
      IF (IDERIV.EQ.0) GO TO 60
      CALL DDERIV (T(IP1), A(IP1MK), KM1, KMIDER, WORK)
C
C *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
C     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      CALL DEVAL (T(IP1), A, KMIDER, WORK, WORK(K+1), WORK(K+K+1), X)
  100 DBUALU = WORK(1)
      RETURN
C
C
  101 CONTINUE
      CALL XERMSG ('SLATEC', 'DBUALU', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  102 CONTINUE
      CALL XERMSG ('SLATEC', 'DBUALU', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DBUALU',
     +   'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
  120 CONTINUE
      CALL XERMSG ('SLATEC', 'DBUALU',
     +   'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
      RETURN
  130 CONTINUE
      CALL XERMSG ('SLATEC', 'DBUALU',
     +   'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
      RETURN
  140 CONTINUE
      CALL XERMSG ('SLATEC', 'DBUALU',
     +   'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
      RETURN
      END
C
      SUBROUTINE DINIT (A, K, AJ)
      INTEGER K, J
      DOUBLE PRECISION A(*), AJ(*)
      DO 30 J=1,K
        AJ(J) = A(J)
   30 CONTINUE
      END
C
      SUBROUTINE DDERIV (T, A, KM1, KMIDER, AJ)
      INTEGER KM1, KMIDER, KMJ, JJ
      DOUBLE PRECISION T(*), A(*), AJ(*), FKMJ
      DO 50 KMJ=KM1,KMIDER,-1
        FKMJ = KMJ
        DO 40 JJ=1,KMJ
          AJ(JJ) = (AJ(JJ+1)-AJ(JJ))/(T(JJ)-T(JJ-KMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
      END
C
      SUBROUTINE DEVAL (T, A, KMIDER, AJ, DP, DM, X)
      INTEGER KMIDER, J, KMJ, ILO
      DOUBLE PRECISION T(*), A(*), AJ(*), DP(*), DM(*), X
      DO 70 J=1,KMIDER
        DP(J) = T(J) - X
        DM(J) = X - T(1-J)
   70 CONTINUE
      DO 90 KMJ=1,KMIDER-1,-1
        ILO = KMJ
        DO 80 JJ=1,KMJ
          AJ(JJ) = (AJ(JJ+1)*DM(ILO)+AJ(JJ)*DP(JJ))/(DM(ILO)+DP(JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
      END

