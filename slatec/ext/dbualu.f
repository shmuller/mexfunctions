*DECK DBUALU
      SUBROUTINE DBUALU (T, A, N, K, D, P, IDERIV, X, INBV, WORK, Y)
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
C          WORK    - work vector of length K*(K+2).
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
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, MFLAG, N,
     1 I1, I2, I3, D, DD, P, PP
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), X, Y(D,P)
C***FIRST EXECUTABLE STATEMENT  DBUALU
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
C
C *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL DINTRV(T, N+1, X, INBV(1), I, MFLAG)
      IF (MFLAG.NE.0) THEN
        IF (MFLAG.EQ.1) THEN
          I = N
        ELSE
          I = K
        END IF
      END IF

   20 IP1 = I + 1
      IP1MK = IP1 - K
      I1 = K + K + 1
      I2 = I1 + K
C
C *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
C     WORK(I) = AJ(I), I=1.K
C     WORK(2*K+I) = T(I)-X, I=1-K.K
C     WORK(3*K+I) = F(I), I=1,K*(K-1)
C
      CALL DINITX (T(IP1), X, K, WORK(I1))
      CALL DINIT3 (WORK(I1), KM1, KMIDER, WORK(I2))
C
      do 40 PP=1,P
        do 30 DD=1,D
          CALL DINIT (A(DD,IP1MK:,PP), K, WORK)
          CALL DEVAL3 (KM1, WORK(I2), WORK)
          Y(DD,PP) = WORK(1)
   30   CONTINUE
   40 CONTINUE
      RETURN
C
C
   99 CONTINUE
      DO 100 I=1,P*D
  100   Y(I,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DINITX (T, X, K, TX)
      INTEGER K, KK
      DOUBLE PRECISION T(*), TX(*), X
      DO 5 KK=1-K,K
        TX(KK) = T(KK) - X
    5 CONTINUE
      END
C
      SUBROUTINE DINIT2 (TX, KM1, KMIDER, F)
      INTEGER KM1, KMIDER, KK, J, I
      DOUBLE PRECISION TX(*), F(*), FACT
      I = 0
      DO 20 KK=KM1,KMIDER,-1
        FACT = KK
        DO 10 J=1,KK
          I = I + 1
          F(I) = FACT/(TX(J)-TX(J-KK))
   10   CONTINUE
   20 CONTINUE
      DO 40 KK=KMIDER-1,1,-1
        DO 30 J=1,KK
          I = I + 1
          F(I) = TX(J)/(TX(J)-TX(J-KK))
   30   CONTINUE
   40 CONTINUE
      END
C
      SUBROUTINE DINIT3 (TX, KM1, KMIDER, F)
      INTEGER KM1, KMIDER, KK, J, I
      DOUBLE PRECISION TX(*), F(*), ONE, FACT, TMP
      ONE = 1
      I = 0
      DO 20 KK=KM1,KMIDER,-1
        FACT = -KK
        DO 10 J=1,KK
          TMP = FACT/(TX(J)-TX(J-KK))
          I = I + 1
          F(I) = TMP
          I = I + 1
          F(I) = -TMP
   10   CONTINUE
   20 CONTINUE
      DO 40 KK=KMIDER-1,1,-1
        DO 30 J=1,KK
          TMP = TX(J)/(TX(J)-TX(J-KK))
          I = I + 1
          F(I) = TMP
          I = I + 1
          F(I) = ONE - TMP
   30   CONTINUE
   40 CONTINUE
      END


      SUBROUTINE DINIT (A, K, AJ)
      INTEGER K, J
      DOUBLE PRECISION A(*), AJ(*)
      DO 5 J=1,K
        AJ(J) = A(J)
    5 CONTINUE
      END
C
      SUBROUTINE DEVAL2 (KM1, KMIDER, F, AJ)
      INTEGER KM1, KMIDER, KK, J, I
      DOUBLE PRECISION F(*), AJ(*), ONE, F1
      ONE = 1
      I = 0
      DO 20 KK=KM1,KMIDER,-1
        DO 10 J=1,KK
          I = I + 1
          AJ(J) = (AJ(J+1)-AJ(J))*F(I)
   10   CONTINUE
   20 CONTINUE
      DO 40 KK=KMIDER-1,1,-1
        DO 30 J=1,KK
          I = I + 1
          F1 = F(I)
          AJ(J) = AJ(J)*F1+AJ(J+1)*(ONE-F1)
   30   CONTINUE
   40 CONTINUE
      END
C
      SUBROUTINE DEVAL3 (KM1, F, AJ)
      INTEGER KM1, KK, J, I
      DOUBLE PRECISION F(*), AJ(*), F1, F2
      I = 0
      DO 20 KK=KM1,1,-1
        DO 10 J=1,KK
          I = I + 1
          F1 = F(I)
          I = I + 1
          F2 = F(I)
          AJ(J) = AJ(J)*F1 + AJ(J+1)*F2
   10   CONTINUE
   20 CONTINUE
      END
      

      SUBROUTINE DDERIV (T, KM1, KMIDER, AJ)
      INTEGER KM1, KMIDER, KMJ, J
      DOUBLE PRECISION T(*), AJ(*), FKMJ
      DO 20 KMJ=KM1,KMIDER,-1
        FKMJ = KMJ
        DO 10 J=1,KMJ
          AJ(J) = (AJ(J+1)-AJ(J))/(T(J)-T(J-KMJ))*FKMJ
   10   CONTINUE
   20 CONTINUE
      END
C
      SUBROUTINE DEVAL (T, X, K, TX, AJ)
      INTEGER K, KK, J, ILO
      DOUBLE PRECISION T(*), TX(*), AJ(*), X
      DO 5 KK=1-K,K
        TX(KK) = T(KK) - X
    5 CONTINUE
      DO 20 KK=K-1,1,-1
        DO 10 J=1,KK
          ILO = J - KK
          AJ(J) = (AJ(J)*TX(J)-AJ(J+1)*TX(ILO))/(TX(J)-TX(ILO))
   10   CONTINUE
   20 CONTINUE
      END

