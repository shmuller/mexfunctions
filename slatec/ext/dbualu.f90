!DECK DBUALU
      SUBROUTINE DBUALU (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, N, &
       I1, I2, I3, D, DD, P, PP, M, MM
      REAL T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
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
  100   Y(I,1,1) = 0.0
      RETURN
      END


      SUBROUTINE DBUALU2 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, &
       WORK2, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, N, &
       I1, I2, I3, D, DD, P, PP, M, MM
      REAL T(*), A(D,N,P), WORK(*), WORK2(D,*), X(*), &
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
  100   Y(I,1,1) = 0.0
      RETURN
      END


      recursive subroutine nd_dot_product(a, sa, b, k, nd, f, res)
      integer sa(nd), k(nd), nd, sa1, sb1, k1, i, j, ndm1
      real a(*), b(*), res(1), f
      k1 = k(1)
      if (nd.gt.1) then
        ndm1 = nd - 1
        sa1 = sa(1)
        sb1 = k1 + 1
        j = 1
        do i=1,k1
          call nd_dot_product(a(j), sa(2), &
                              b(sb1), k(2), ndm1, f*b(i), res)
          j = j + sa1
        end do
      else
        res(1) = res(1) + f*dot_product(a(1:k1), b(1:k1))
      end if
      end

      recursive subroutine nd_dot_product2(a, sa, b, sb, k, nd, f, res)
      integer sa(nd), sb(nd), k(nd), nd, sa1, sb1, k1, i, j, ndm1
      real a(*), b(*), res(1), f
      k1 = k(1)
      if (nd.gt.1) then
        ndm1 = nd - 1
        sa1 = sa(1)
        sb1 = sb(1)
        j = 1
        do i=1,k1
          call nd_dot_product2(a(j), sa(2), &
                               b(sb1), sb(2), k(2), ndm1, f*b(i), res)
          j = j + sa1
        end do
      else
        res(1) = res(1) + f*dot_product(a(1:k1), b(1:k1))
      end if
      end


      SUBROUTINE DBDER(T, K, IDERIV, X, VNIKX)
      INTEGER K, IDERIV, KMIDER, L, J
      REAL T(*), X, VNIKX(*), VM, VMPREV, FKMJ, DL, DR
      KMIDER = K - IDERIV
      VNIKX(1) = 1.0
      DO 20 J=1,KMIDER-1
        VMPREV = 0.0
        DO 10 L=1,J
          DR = T(L) - X
          DL = X - T(L-J)
          VM = VNIKX(L)/(DR+DL)
          VNIKX(L) = VMPREV + VM*DR
          VMPREV = VM*DL
   10   CONTINUE
        VNIKX(J+1) = VMPREV
   20 CONTINUE
      DO 40 J=KMIDER,K-1
        FKMJ = J
        VMPREV = 0.0
        DO 30 L=1,J
          VM = VNIKX(L)/(T(L)-T(L-J))*FKMJ
          VNIKX(L) = VMPREV - VM
          VMPREV = VM
   30   CONTINUE
        VNIKX(J+1) = VMPREV
   40 CONTINUE
      END


      SUBROUTINE DBSPVN2 (T, K, X, VNIKX, WORK)
      INTEGER JP1, JP1ML, K, J, L
      REAL T(*), VM, VMPREV, VNIKX(*), WORK(*), X
!     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
!     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
!***FIRST EXECUTABLE STATEMENT  DBSPVN2
      VNIKX(1) = 1.0
      DO 40 J=1,K-1
        WORK(J) = T(J) - X
        WORK(K+J) = X - T(1-J)
        VMPREV = 0.0
        JP1 = J + 1
        DO 30 L=1,J
          JP1ML = JP1 - L
          VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
          VNIKX(L) = VM*WORK(L) + VMPREV
          VMPREV = VM*WORK(K+JP1ML)
   30   CONTINUE
        VNIKX(JP1) = VMPREV
   40 CONTINUE
      END


      SUBROUTINE DBUALND (NDIM, T, A, N, K, S, IDERIV, X, M, INBV, &
       WORK, Y)
      INTEGER NDIM, N(NDIM), K(NDIM), S(NDIM), IDERIV(NDIM), &
       INBV(NDIM), I, OFFS, D, ND, KD, JT, JB, M, MM
      REAL T(*), A(*), WORK(*), X(NDIM,M), Y(M), F
!***FIRST EXECUTABLE STATEMENT  DBUAL
      if (ANY(K.LE.IDERIV)) GO TO 99

      S(NDIM) = 1
      DO 10 D=NDIM,2,-1
        S(D-1) = N(D)*S(D)
   10 CONTINUE

      DO 50 MM=1,M
        JT = 1
        JB = 1
        OFFS = 1
        DO 20 D=1,NDIM
          ND = N(D)
          KD = K(D)
          IF (MM.GT.1 .AND. X(D,MM-1).EQ.X(D,MM)) THEN
            I = INBV(D)
          ELSE
            CALL DFINDI (T(JT), ND, KD, X(D,MM), INBV(D), I)
            CALL DBDER (T(JT+I), KD, IDERIV(D), X(D,MM), WORK(JB))
          END IF
          JT = JT + ND + KD
          JB = JB + KD
          OFFS = OFFS + (I-KD)*S(D)
   20   CONTINUE

        F = 1
        Y(MM) = 0
        CALL nd_dot_product(A(OFFS), S, WORK, K, NDIM, F, Y(MM))
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 MM=1,M
  100   Y(MM) = 0.0
      RETURN
      END


      SUBROUTINE DBSPGD (NDIM, T, N, K, S, IDERIV, X, M, I, B)
      INTEGER NDIM, N(*), K(*), S(*), IDERIV(*), M(*), I(*), &
       JT, JI, JB, J, D, ND, KD, SD, INBV, IDER
      REAL T(*), X(*), B(*)
      S(NDIM) = 1
      DO 5 D=NDIM,2,-1
        S(D-1) = N(D)*S(D)
    5 CONTINUE
      JT = 1
      JI = 1
      JB = 1
      DO 15 D=1,NDIM
        ND = N(D)
        KD = K(D)
        SD = S(D)
        IDER = IDERIV(D)
        INBV = 1
        DO 10 J=1,M(D)
          CALL DFINDI (T(JT), ND, KD, X(JI), INBV, I(JI))
          CALL DBDER (T(JT+I(JI)), KD, IDER, X(JI), B(JB))
          I(JI) = SD*(I(JI)-KD)
          JI = JI + 1
          JB = JB + KD
   10   CONTINUE
        JT = JT + ND + KD
   15 CONTINUE
      END

      RECURSIVE SUBROUTINE LOOPGD (A, SA, I, B, SB, SSB, DSB, K, M, &
       ND, D, R, IR)
      INTEGER I(*), SA(*), SB(*), SSB(*), DSB(*), K(*), M(*), ND, D, &
       IR, MD, KD, J
      REAL A(*), B(*), R(*), F
      IF (D.LE.ND) THEN
        MD = M(D)
        KD = K(D)
        SSB(D) = SB(D)
        DO 20 J=1,MD
          CALL LOOPGD (A(I(J)+1), SA, I(MD+1), B, SB, SSB, &
                       DSB, K, M, ND, D+1, R, IR)
          SSB(D) = SSB(D)+KD
   20   CONTINUE
      ELSE
        DO 25 J=1,ND-1
          DSB(J) = SSB(J+1) - SSB(J) + 1
   25   CONTINUE
        F = 1.0
        R(IR) = 0.0
        CALL nd_dot_product2(A, SA, B(SSB(1)+1), DSB, K, ND, F, R(IR))
        IR = IR + 1
      END IF
      END


      SUBROUTINE DBUALGD (NDIM, T, A, N, K, S, IDERIV, X, M, I, B, R)
      INTEGER NDIM, N(NDIM), K(NDIM), S(4*NDIM), IDERIV(NDIM), &
       M(NDIM), D, I(*), IR
      REAL T(*), A(*), X(*), B(*), R(*)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      CALL DBSPGD (NDIM, T, N, K, S, IDERIV, X, M, I, B)
      S(NDIM+1) = 0
      DO 30 D=1,NDIM-1
        S(NDIM+D+1) = S(NDIM+D) + K(D)*M(D)
   30 CONTINUE
      IR = 1
      CALL LOOPGD (A, S, I, B, S(NDIM+1), S(2*NDIM+1), S(3*NDIM+1), &
                   K, M, NDIM, 1, R, IR)
      END


      SUBROUTINE DBUAL3D (T, A, N, K, IDERIV, X, M, I, B, R)
      INTEGER N(3), K(3), S(3), IDERIV(3), M(3), D, I(*), IR, IX, IXY, &
       JX, JY, JZ, LX, LY, LZ, SBX, SBY, SBZ, I1, I12, &
       KX, KY, KZ, NZ, NYZ, MX, MXY, MXYZ, MKX, MKXY
      REAL T(*), A(*), X(*), B(*), R(*), BX, BXY, RI
!***FIRST EXECUTABLE STATEMENT  DBUAL
      CALL DBSPGD (3, T, N, K, S, IDERIV, X, M, I, B)
      KX = K(1)
      KY = K(2)
      KZ = K(3)
      NZ = N(3)
      NYZ = N(2)*NZ
      MX = M(1)
      MXY = MX + M(2)
      MXYZ = MXY + M(3)
      MKX = MX*KX
      MKXY = MKX + M(2)*KY

      IR = 1
      SBX = 0
      DO JX=1,MX
        IX = I(JX)
        SBY = MKX
        DO JY=MX+1,MXY
          IXY = IX + I(JY)
          SBZ = MKXY
          DO JZ=MXY+1,MXYZ

            RI = 0.0D0
            I1 = IXY + I(JZ)
            DO LX=1,KX
              I12 = I1
              BX = B(SBX+LX)
              DO LY=1,KY
                BXY = BX*B(SBY+LY)
                DO LZ=1,KZ
                  RI = RI + A(I12+LZ)*BXY*B(SBZ+LZ)
                END DO
                I12 = I12 + NZ
              END DO
              I1 = I1 + NYZ
            END DO
            R(IR) = RI
            IR = IR + 1

            SBZ = SBZ + KZ
          END DO
          SBY = SBY + KY
        END DO
        SBX = SBX + KX
      END DO
      END


      SUBROUTINE DBUAL (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, N, &
       D, P, PP, M, MM
      REAL T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      IF (K.LE.IDERIV) GO TO 99
      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K
        CALL DBDER (T(IP1), K, IDERIV, X(MM), WORK)
        DO 20 PP=1,P
           Y(:,MM,PP) = MATMUL(A(:,IP1MK:I,PP), WORK(1:K))
   20   CONTINUE
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0
      RETURN
      END


      SUBROUTINE DBUAL2 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), K, IMK, J, N, &
       D, P, PP, M, MM
      REAL T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      IF (K.LE.IDERIV) GO TO 99
      Y = 0.0
      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IMK = I - K
        CALL DBDER (T(I+1), K, IDERIV, X(MM), WORK)
        DO 20 PP=1,P
          DO 10 J=1,K
            Y(:,MM,PP) = Y(:,MM,PP) + A(:,IMK+J,PP)*WORK(J)
   10     CONTINUE
   20   CONTINUE
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0
      RETURN
      END


      SUBROUTINE DBUAL3 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, &
       WORK2, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, N, D, P, M, MM
      REAL T(*), A(D,N,P), WORK(*), WORK2(D,P,*), X(*), &
       Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      IF (K.LE.IDERIV) GO TO 99
      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K
        CALL DBDER (T(IP1), K, IDERIV, X(MM), WORK)

        WORK2(:,:,1:K) = RESHAPE( &
          A(:,IP1MK:I,:), (/D,P,K/), ORDER=(/1,3,2/))

        Y(:,MM,:) = RESHAPE(MATMUL(RESHAPE( &
          WORK2(:,:,1:K), (/D*P,K/)), WORK(1:K)), (/D,P/))
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,P*M*D
  100   Y(I,1,1) = 0.0
      RETURN
      END


      SUBROUTINE DBUAL4 (T, A, N, K, D, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, N, D, M, MM
      REAL T(*), A(D,N), WORK(*), X(*), Y(D,M)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      IF (K.LE.IDERIV) GO TO 99
      DO 50 MM=1,M
        CALL DFINDI (T, N, K, X(MM), INBV(1), I)
        IP1 = I + 1
        IP1MK = IP1 - K
        CALL DBDER (T(IP1), K, IDERIV, X(MM), WORK)
        Y(:,MM) = MATMUL(A(:,IP1MK:I), WORK(1:K))
   50 CONTINUE
      RETURN

   99 CONTINUE
      DO 100 I=1,M*D
  100   Y(I,1) = 0.0
      RETURN
      END


      SUBROUTINE DBVAL1 (T, A, N, K, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, KMIDER, KM1, N, &
       I1, I2, I3, M, MM
      REAL T(*), A(N), WORK(*), X(*), Y(M)
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
  100   Y(I) = 0.0
      RETURN
      END


      SUBROUTINE DBVALI (T, A, N, K, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), K, KMIDER, KM1, N, NP1, &
       M, MM, MFLAG, KK, J, IMK, IPJ
      REAL T(*), A(N), WORK(*), X(*), Y(M), F1, F2, FKMJ
!***FIRST EXECUTABLE STATEMENT  DBVAL1
      KMIDER = K - IDERIV
      IF (KMIDER.LE.0) GO TO 99
      KM1 = K - 1
      NP1 = N + 1
      DO 50 MM=1,M
        CALL INTRV(T, NP1, X(MM), INBV(1), I, MFLAG)
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
  100   Y(I) = 0.0
      RETURN
      END


      SUBROUTINE DFINDI(T, N, K, X, INBV, I)
! *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
!     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      INTEGER N, K, INBV, I, MFLAG
      REAL T(*), X
      CALL INTRV(T, N+1, X, INBV, I, MFLAG)
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
      REAL T(*), TX(*), X
      DO 5 KK=1-K,K
        TX(KK) = T(KK) - X
    5 CONTINUE
      END

      SUBROUTINE DINIT2 (TX, KM1, KMIDER, F)
      INTEGER KM1, KMIDER, KK, J, I
      REAL TX(*), F(*), FACT
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
      REAL TX(*), F(*), ONE, FACT, TMP
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
      REAL AJ(*), A(D,*)
      DO 5 J=1,K
        AJ(J) = A(DD,J)
    5 CONTINUE
      END

      SUBROUTINE DEVAL2 (KM1, KMIDER, F, AJ)
      INTEGER KM1, KMIDER, KK, J, I
      REAL F(*), AJ(*), ONE, F1
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
      REAL F(*), AJ(*), F1, F2
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
      REAL F(*), AJ(D,*), F1, F2
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
      REAL T(*), AJ(*), FKMJ
      DO 20 KMJ=KM1,KMIDER,-1
        FKMJ = KMJ
        DO 10 J=1,KMJ
          AJ(J) = (AJ(J+1)-AJ(J))/(T(J)-T(J-KMJ))*FKMJ
   10   CONTINUE
   20 CONTINUE
      END

      SUBROUTINE DEVAL (T, X, K, TX, AJ)
      INTEGER K, KK, J, ILO
      REAL T(*), TX(*), AJ(*), X
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


      SUBROUTINE APTKNT (X, N, K, T)
      INTEGER N, K, I, KM1
      REAL X(N), T(N+K), FKM1
      KM1 = K-1
      FKM1 = KM1
      DO 10 I=1,K
        T(I) = X(1)
        T(N+I) = X(N)
   10 CONTINUE
      DO 20 I=K+1,N
        T(I) = SUM(X(I-KM1:I-1))/FKM1
   20 CONTINUE
      END
