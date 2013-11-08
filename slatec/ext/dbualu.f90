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


      recursive subroutine nd_dot_product(a, s, b, n, nd, f, res)
      integer s(nd), n(nd), nd, s1, n1, i, j, n1p1, ndm1
      double precision a(*), b(*), res(1), f
      s1 = s(1)
      n1 = n(1)
      n1p1 = n1 + 1
      ndm1 = nd - 1
      if (nd.gt.1) then
        j = 1
        do i=1,n1
          call nd_dot_product(a(j), s(2), &
                              b(n1p1), n(2), ndm1, f*b(i), res)
          j = j + s1
        end do
      else
        res(1) = res(1) + f*dot_product(a(1:n1), b(1:n1))
      end if
      end


      SUBROUTINE DBDER(T, K, IDERIV, X, VNIKX)
      INTEGER K, IDERIV, KMIDER, L, J
      DOUBLE PRECISION T(*), X, VNIKX(*), VM, VMPREV, FKMJ, DL, DR
      KMIDER = K - IDERIV
      VNIKX(1) = 1.0D0
      DO 20 J=1,KMIDER-1
        VMPREV = 0.0D0
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
        VMPREV = 0.0D0
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
      DOUBLE PRECISION T(*), VM, VMPREV, VNIKX(*), WORK(*), X
!     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
!     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
!***FIRST EXECUTABLE STATEMENT  DBSPVN2
      VNIKX(1) = 1.0D0
      DO 40 J=1,K-1
        WORK(J) = T(J) - X
        WORK(K+J) = X - T(1-J)
        VMPREV = 0.0D0
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
      DOUBLE PRECISION T(*), A(*), WORK(*), X(NDIM,M), Y(M), F
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
  100   Y(MM) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUALGD (NDIM, T, A, N, K, SI, SB, SA, SR, &
       X, M, I, B, R)
      INTEGER NDIM, N(NDIM), K(NDIM), SI(NDIM), SB(NDIM), &
       SA(NDIM), SR(NDIM), M(NDIM), ND, KD, &
       D, J, I(*), JX, JY, JZ, INBV(3), JT, JB, JI, &
       IX, IY, IZ, LX, LY, LZ
      DOUBLE PRECISION T(*), A(*), X(*), XJ, BX, BY, BZ, &
       B(*), R(*), S
!***FIRST EXECUTABLE STATEMENT  DBUAL
      SI(1) = 0
      SB(1) = 0
      DO 5 D=1,NDIM-1
        SI(D+1) = SI(D) + M(D)
        SB(D+1) = SB(D) + K(D)*M(D)
    5 CONTINUE
      INBV = 1
      JT = 1
      JI = 1
      JB = 1
      DO 15 D=1,NDIM
        ND = N(D)
        KD = K(D)
        DO 10 J=1,M(D)
          CALL DFINDI (T(JT), ND, KD, X(JI), INBV(D), I(JI))
          CALL DBDER (T(JT+I(JI)), KD, 0, X(JI), B(JB))
          JI = JI + 1
          JB = JB + KD
   10   CONTINUE
        JT = JT + ND + KD
   15 CONTINUE

      SA(NDIM) = 1
      SR(NDIM) = 1
      DO 16 D=NDIM,2,-1
        SA(D-1) = N(D)*SA(D)
        SR(D-1) = M(D)*SR(D)
   16 CONTINUE
      DO 22 JX=1,M(1)
        DO 21 JY=1,M(2)
          DO 20 JZ=1,M(3)
            S = 0.0D0
            DO 19 LX=1,K(1)
              DO 18 LY=1,K(2)
                DO 17 LZ=1,K(3)
                  IX = I(SI(1)+JX) - K(1)
                  IY = I(SI(2)+JY) - K(2)
                  IZ = I(SI(3)+JZ) - K(3)
                  BX = B(SB(1)+(JX-1)*K(1)+LX)
                  BY = B(SB(2)+(JY-1)*K(2)+LY)
                  BZ = B(SB(3)+(JZ-1)*K(3)+LZ)
                  S = S &
                    + A(SA(1)*(IX+LX-1)+SA(2)*(IY+LY-1)+IZ+LZ)*BX*BY*BZ
   17           CONTINUE
   18         CONTINUE
   19       CONTINUE
            R(SR(1)*(JX-1)+SR(2)*(JY-1)+JZ) = S
   20     CONTINUE
   21   CONTINUE
   22 CONTINUE
      END


      SUBROUTINE DBUAL3D (NDIM, T, A, N, K, SI, SB, NX, NY, NZ, &
       X, M, MX, MY, MZ, I, B, AX, AXY, R)
      INTEGER NDIM, N(NDIM), K(NDIM), SI(NDIM), SB(NDIM), M(NDIM), &
       NX, NY, NZ, KX, KY, KZ, ND, KD, &
       D, J, MX, MY, MZ, I(*), JX, JY, JZ, INBV(3), JT, JB, JI, &
       IX, IY, IZ
      DOUBLE PRECISION T(*), A(NZ,NY,NX), X(*), XJ, &
       B(*), AX(NZ,NY), AXY(NZ), R(MZ,MY,MX)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      SI(1) = 0
      SB(1) = 0
      DO 5 D=1,NDIM-1
        SI(D+1) = SI(D) + M(D)
        SB(D+1) = SB(D) + K(D)*M(D)
    5 CONTINUE
      INBV = 1
      JT = 1
      JI = 1
      JB = 1
      DO 15 D=1,NDIM
        ND = N(D)
        KD = K(D)
        DO 10 J=1,M(D)
          CALL DFINDI (T(JT), ND, KD, X(JI), INBV(D), I(JI))
          CALL DBDER (T(JT+I(JI)), KD, 0, X(JI), B(JB))
          JI = JI + 1
          JB = JB + KD
   10   CONTINUE
        JT = JT + ND + KD
   15 CONTINUE
      KX = K(1)
      KY = K(2)
      KZ = K(3)
      DO 22 JX=1,MX
        IX = I(SI(1)+JX)
        AX(:,:) = RESHAPE(MATMUL(RESHAPE(A(:,:,IX+1-KX:IX), &
                  (/NZ*NY,KX/)), &
                  B(SB(1)+(JX-1)*KX+1:SB(1)+JX*KX)), (/NZ,NY/))
        DO 21 JY=1,MY
          IY = I(SI(2)+JY)
          AXY(:) = MATMUL(AX(:,IY+1-KY:IY), &
                          B(SB(2)+(JY-1)*KY+1:SB(2)+JY*KY))
          DO 20 JZ=1,MZ
            IZ = I(SI(3)+JZ)
            R(JZ,JY,JX) = DOT_PRODUCT(AXY(IZ+1-KZ:IZ), &
                                      B(SB(3)+(JZ-1)*KZ+1:SB(3)+JZ*KZ))
   20     CONTINUE
   21   CONTINUE
   22 CONTINUE
      END


      SUBROUTINE DBUAL (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, N, &
       D, P, PP, M, MM
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
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
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUAL2 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), K, IMK, J, N, &
       D, P, PP, M, MM
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), X(*), Y(D,M,P)
!***FIRST EXECUTABLE STATEMENT  DBUAL
      IF (K.LE.IDERIV) GO TO 99
      Y = 0.0D0
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
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUAL3 (T, A, N, K, D, P, IDERIV, X, M, INBV, WORK, &
       WORK2, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, N, D, P, M, MM
      DOUBLE PRECISION T(*), A(D,N,P), WORK(*), WORK2(D,P,*), X(*), &
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
  100   Y(I,1,1) = 0.0D0
      RETURN
      END


      SUBROUTINE DBUAL4 (T, A, N, K, D, IDERIV, X, M, INBV, WORK, Y)
      INTEGER I, IDERIV, INBV(*), IP1, K, IP1MK, N, D, M, MM
      DOUBLE PRECISION T(*), A(D,N), WORK(*), X(*), Y(D,M)
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

