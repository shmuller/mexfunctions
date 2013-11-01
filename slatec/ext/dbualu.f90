!DECK DBUALU
      SUBROUTINE DBUALU (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, &
       WORK2, Y)
!***BEGIN PROLOGUE  DBUALU
!***PURPOSE  Evaluate the B-representation of a B-spline at X for the
!            function value or any of its derivatives.
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      DOUBLE PRECISION (BVALU-S, DBUALU-D)
!***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract   **** a double precision routine ****
!         DBUALU is the BVALUE function of the reference.
!
!         DBUALU evaluates the B-representation (T,A,N,K) of a B-spline
!         at X for the function value on IDERIV=0 or any of its
!         derivatives on IDERIV=1,2,...,K-1.  Right limiting values
!         (right derivatives) are returned except at the right end
!         point X=T(N+1) where left limiting values are computed.  The
!         spline is defined on T(K) .LE. X .LE. T(N+1).  DBUALU returns
!         a fatal error message when X is outside of this interval.
!
!         To compute left derivatives or left limiting values at a
!         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!         DBUALU calls DINTRV
!
!     Description of Arguments
!
!         Input      T,A,X are double precision
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K .GE. 1
!          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!                    IDERIV = 0 returns the B-spline value
!          X       - argument, T(K) .LE. X .LE. T(N+1)
!          INBV    - an initialization parameter which must be set
!                    to 1 the first time DBUALU is called.
!
!         Output     WORK,DBUALU are double precision
!          INBV    - INBV contains information for efficient process-
!                    ing after the initial call and INBV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INBV parameters.
!          WORK    - work vector of length K*(K+2), used as:
!                      WORK(I) = AJ(I), I=1.K
!                      WORK(2*K+I) = T(I)-X, I=1-K.K
!                      WORK(3*K+I) = F(I), I=1,K*(K-1)
!          DBUALU  - value of the IDERIV-th derivative at X
!
!     Error Conditions
!         An improper input is a fatal error
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  DINTRV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBUALU
!
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, N, &
       I1, I2, I3, D, DD, P, PP, M, MM
      DOUBLE PRECISION T(*), A(D,N,*), WORK(*), WORK2(D,*), X(*), &
       Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUALU
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99

      KM1 = K - 1
      I1 = 1
      I2 = I1 + K
      I3 = I2 + K

      do 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K

!        CALL DINITX (T(IP1), X(MM), K, WORK(I2))
        WORK(I1:I3-1) = T(IP1MK:I+K) - X(MM)

        CALL DINIT3 (WORK(I2), KM1, KMIDER, WORK(I3))

        DO 40 PP=1,P
          WORK2(:,1:K) = A(:,IP1MK:I,PP)
          CALL DEVAL3 (KM1, WORK(I3), WORK2, D)
          Y(:,MM,PP) = WORK2(:,1)
   40   CONTINUE

!        do 40 PP=1,P
!          do 30 DD=1,D
!!            CALL DINIT (A(:,IP1MK:,PP), K, WORK, D, DD)
!            WORK(1:K) = A(DD,IP1MK:I,PP)
!            CALL DEVAL3 (KM1, WORK(I3), WORK)
!            Y(DD,MM,PP) = WORK(1)
!   30     CONTINUE
!   40   CONTINUE
   50 CONTINUE
      RETURN


   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0D0
      RETURN
      END

      SUBROUTINE DFINDI(T, N, K, X, INBV, I)
! *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
!     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      INTEGER N, K, INBV, I, MFLAG
      DOUBLE PRECISION T(*), X
      CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (MFLAG.NE.0) THEN
        IF (MFLAG.EQ.1) THEN
          I = N
        ELSE
          I = K
        END IF
      END IF
      END

      SUBROUTINE DINITX (T, X, K, TX)
      INTEGER K, KK
      DOUBLE PRECISION T(*), TX(*), X
      DO 5 KK=1-K,K
        TX(KK) = T(KK) - X
    5 CONTINUE
      END

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


      SUBROUTINE DINIT (A, K, AJ, D, DD)
      INTEGER K, J, D, DD
      DOUBLE PRECISION AJ(*), A(D,*)
      DO 5 J=1,K
        AJ(J) = A(DD,J)
    5 CONTINUE
      END

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

      SUBROUTINE DEVAL3 (KM1, F, AJ, D)
      INTEGER KM1, KK, J, I, D
      DOUBLE PRECISION F(*), AJ(D,*), F1, F2
      I = 0
      DO 20 KK=KM1,1,-1
        DO 10 J=1,KK
          I = I + 1
          F1 = F(I)
          I = I + 1
          F2 = F(I)
          AJ(:,J) = AJ(:,J)*F1 + AJ(:,J+1)*F2
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

