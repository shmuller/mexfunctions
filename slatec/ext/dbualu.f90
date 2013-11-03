!DECK DBUALU
      SUBROUTINE DBUALU (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, N, &
       I1, I2, I3, D, DD, P, PP, M, MM
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUALU
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99

      KM1 = K - 1
      I1 = K + 1
      I2 = I1 + K
      I3 = I2 + K

      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K

        CALL DINITX (T(IP1), X(MM), K, WORK(I2))
        WORK(I1:I3-1) = T(IP1MK:I+K) - X(MM)

        CALL DINIT3 (WORK(I2), KM1, KMIDER, WORK(I3))

        do 40 PP=1,P
          do 30 DD=1,D
            WORK(1:K) = A(DD,IP1MK:I,PP)
            CALL DEVAL3 (KM1, WORK(I3), WORK)
            Y(DD,MM,PP) = WORK(1)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUALU2 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, &
       WORK2, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, N, &
       I1, I2, I3, D, DD, P, PP, M, MM
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), WORK2(D,*), X(*), &
       Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUALU
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99

      KM1 = K - 1
      I1 = 1
      I2 = I1 + K
      I3 = I2 + K

      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K

        CALL DINITX (T(IP1), X(MM), K, WORK(I2))
        WORK(I1:I3-1) = T(IP1MK:I+K) - X(MM)

        CALL DINIT3 (WORK(I2), KM1, KMIDER, WORK(I3))

        DO 40 PP=1,P
          WORK2(:,1:K) = A(:,IP1MK:I,PP)
          CALL DEVAL3V (KM1, WORK(I3), WORK2, D)
          Y(:,MM,PP) = WORK2(:,1)
   40   CONTINUE
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUAL (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KP1, N, &
       D, P, PP, M, MM, IWORK
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
      KP1 = K + 1

      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K

        CALL DBSPVN (T, K, K, 1, X(MM), I, WORK, WORK(KP1), IWORK)
        DO 20 PP=1,P
           Y(:,MM,PP) = MATMUL(A(:,IP1MK:I,PP), WORK(1:K))
   20   CONTINUE
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUAL2 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IMK, J, KMIDER, KP1, N, &
       D, P, PP, M, MM, IWORK
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
      KP1 = K + 1

      Y = 0.0D0
      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IMK = I - K

        CALL DBSPVN (T, K, K, 1, X(MM), I, WORK, WORK(KP1), IWORK)
        DO 20 PP=1,P
          DO 10 J=1,K
            Y(:,MM,PP) = Y(:,MM,PP) + A(:,IMK+J,PP)*WORK(J)
   10     CONTINUE
   20   CONTINUE
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUAL3 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, &
       WORK2, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, J, KMIDER, KP1, N, &
       D, P, PP, M, MM, IWORK
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), WORK2(D,P,*), X(*), &
       Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
      KP1 = K + 1

      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K

        CALL DBSPVN (T, K, K, 1, X(MM), I, WORK, WORK(KP1), IWORK)

        WORK2(:,:,1:K) = RESHAPE( &
          A(:,IP1MK:I,:), (/D,P,K/), ORDER=(/1,3,2/))

        Y(:,MM,:) = RESHAPE(MATMUL(RESHAPE( &
          WORK2(:,:,1:K), (/D*P,K/)), WORK(1:K)), (/D,P/))
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUAL4 (T, A, N, K, D, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KP1, N, &
       D, P, PP, M, MM, IWORK
      DOUBLE PRECISION T(*), A(D,N), WORK(*), X(*), Y(D,M)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
      KP1 = K + 1

      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K

        CALL DBSPVN (T, K, K, 1, X(MM), I, WORK, WORK(KP1), IWORK)
        Y(:,MM) = MATMUL(A(:,IP1MK:I), WORK(1:K))
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,M*D
  100   Y(I,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBVAL1 (T, A, N, K, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, N, &
       I1, I2, I3, M, MM
      DOUBLE PRECISION T(*), A(N), WORK(*), X(*), Y(M)
!***FIRST EXECUTABLE STATEMENT  DBVAL1
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99

      KM1 = K - 1
      I1 = K + 1
      I2 = I1 + K
      I3 = I2 + K

      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K

        WORK(1:K) = A(IP1MK:I)
        CALL DDERIV (T(IP1), KM1, KMIDER, WORK)
        CALL DEVAL (T(IP1), X(MM), KMIDER, WORK(I2), WORK)
        Y(MM) = WORK(1)
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,M
  100   Y(I) = 0.0D0
      RETURN
      END


      SUBROUTINE DBVALI (T, A, N, K, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), K, KMIDER, KM1, N, NP1, &
       M, MM, MFLAG, KK, J, IMK, IPJ
      DOUBLE PRECISION T(*), A(N), WORK(*), X(*), Y(M), F1, F2, FKMJ
!***FIRST EXECUTABLE STATEMENT  DBVAL1
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
      KM1 = K - 1
      NP1 = N + 1
      DO 50 MM=1,M
        CALL DINTRV(T, NP1, X(MM), INBV(1), I, MFLAG)
        IF (MFLAG.NE.0) THEN
          IF (MFLAG.EQ.1) THEN
            I = N
          ELSE
            I = K
          END IF
        END IF
        IMK = I - K
        DO 5 KK=1,K
          WORK(KK) = A(IMK+KK)
    5   CONTINUE
        DO 20 KK=KM1,KMIDER,-1
          FKMJ = KK
          DO 10 J=1,KK
            IPJ = I + J
            WORK(J) = (WORK(J+1)-WORK(J))/(T(IPJ)-T(IPJ-KK))*FKMJ
   10     CONTINUE
   20   CONTINUE
        DO 40 KK=KMIDER-1,1,-1
          DO 30 J=1,KK
            IPJ = I + J
            F1 = T(IPJ) - X(MM)
            F2 = T(IPJ-KK) - X(MM)
            WORK(J) = (WORK(J)*F1-WORK(J+1)*F2)/(F1-F2)
   30     CONTINUE
   40   CONTINUE
        Y(MM) = WORK(1)
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,M
  100   Y(I) = 0.0D0
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

      SUBROUTINE DEVAL3V (KM1, F, AJ, D)
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

