*DECK CGBQC
      SUBROUTINE CGBQC (LUN, KPRINT, NERR)
C***BEGIN PROLOGUE  CGBQC
C***PURPOSE  Quick check for CGBFA, CGBCO, CGBSL and CGBDI.
C***LIBRARY   SLATEC
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C    LET  A*X=B  BE A COMPLEX LINEAR SYSTEM WHERE THE MATRIX  A  IS
C    OF THE PROPER TYPE FOR THE LINPACK SUBROUTINES BEING TESTED.
C    THE VALUES OF  A  AND  B  AND THE PRE-COMPUTED VALUES OF  C
C    (THE SOLUTION VECTOR),  DC  (DETERMINANT OF  A ), AND
C    RCND  (RCOND) ARE ENTERED WITH DATA STATEMENTS.
C
C    THE COMPUTED TEST RESULTS FOR  X,  RCOND  AND THE DETER-
C    MINANT ARE COMPARED TO THE STORED PRE-COMPUTED VALUES.
C    FAILURE OF THE TEST OCCURS WHEN AGREEMENT TO 3 SIGNIFICANT
C    DIGITS IS NOT ACHIEVED AND AN ERROR MESSAGE INDICATING WHICH
C    LINPACK SUBROUTINE FAILED AND WHICH QUANTITY WAS INVOLVED IS
C    PRINTED.  A SUMMARY LINE IS ALWAYS PRINTED.
C
C    NO INPUT ARGUMENTS ARE REQUIRED.
C    ON RETURN,  NERR  (INTEGER TYPE) CONTAINS THE TOTAL COUNT OF
C    ALL FAILURES DETECTED BY CGBQC.
C
C***ROUTINES CALLED  CGBCO, CGBDI, CGBFA, CGBSL
C***REVISION HISTORY  (YYMMDD)
C   801015  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Restructured using IF-THEN-ELSE-ENDIF, moved an ARITHMETIC
C           STATEMENT FUNCTION ahead of the FIRST EXECUTABLE STATEMENT
C           record and cleaned up FORMATs.  (RWC)
C***END PROLOGUE  CGBQC
      COMPLEX ABD(6,4),AT(7,4),B(4),BT(4),C(4),DET(2),DC(2),
     1 Z(4),XA,XB
      REAL R,RCOND,RCND,DELX
      CHARACTER KFAIL*39,KPROG*19
      INTEGER LDA,N,IPVT(4),INFO,I,J,INDX,NERR
      INTEGER ML,MU
      DATA ABD/(0.E0,0.E0),(0.E0,0.E0),(0.E0,0.E0),(0.E0,0.E0),
     1 (2.E0,0.E0),(0.E0,1.E0),(0.E0,0.E0),(0.E0,0.E0),
     2 (0.E0,0.E0),(0.E0,-1.E0),(2.E0,0.E0),(0.E0,0.E0),
     3 (0.E0,0.E0),(0.E0,0.E0),(0.E0,0.E0),(0.E0,0.E0),
     4 (3.E0,0.E0),(0.E0,1.E0),(0.E0,0.E0),(0.E0,0.E0),
     5 (0.E0,0.E0),(0.E0,-1.E0),(4.E0,0.E0),(0.E0,0.E0)/
      DATA B/(3.E0,2.E0),(-1.E0,3.E0),(0.E0,-4.E0),(5.E0,0.E0)/
      DATA C/(1.E0,1.E0),(0.E0,1.E0),(0.E0,-1.E0),(1.E0,0.E0)/
      DATA DC/(3.3E0,0.E0),(1.0E0,0.E0)/
      DATA KPROG/'GBFA GBCO GBSL GBDI'/
      DATA KFAIL/'INFO RCOND SOLUTION DETERMINANT INVERSE'/
      DATA RCND/.24099E0/
C
      DELX(XA,XB)=ABS(REAL(XA-XB))+ABS(AIMAG(XA-XB))
C***FIRST EXECUTABLE STATEMENT  CGBQC
      LDA = 7
      N = 4
      ML = 1
      MU = 3
      NERR = 0
C
C     FORM AT FOR CGBFA AND BT FOR CGBSL, TEST CGBFA
C
      DO 20 J=1,N
         BT(J) = B(J)
         DO 10 I=1,6
            AT(I,J) = ABD(I,J)
   10    CONTINUE
   20 CONTINUE
C
      CALL CGBFA(AT,LDA,N,ML,MU,IPVT,INFO)
      IF (INFO .NE. 0) THEN
         WRITE (LUN,201) KPROG(1:4),KFAIL(1:4)
         NERR = NERR + 1
      ENDIF
C
C     TEST CGBSL FOR JOB=0
C
      CALL CGBSL(AT,LDA,N,ML,MU,IPVT,BT,0)
      INDX = 0
      DO 40 I=1,N
         IF (DELX(C(I),BT(I)) .GT. .0001) INDX=INDX+1
   40 CONTINUE
C
      IF (INDX .NE. 0) THEN
         WRITE (LUN,201) KPROG(11:14),KFAIL(12:19)
         NERR = NERR + 1
      ENDIF
C
C     FORM AT FOR CGBCO AND BT FOR CGBSL, TEST CGBCO
C
      DO 70 J=1,N
         BT(J) = B(J)
         DO 60 I=1,6
            AT(I,J) = ABD(I,J)
   60    CONTINUE
   70 CONTINUE
C
      CALL CGBCO(AT,LDA,N,ML,MU,IPVT,RCOND,Z)
      R = ABS(RCND-RCOND)
      IF (R .GE. .0001) THEN
         WRITE (LUN,201) KPROG(6:9),KFAIL(6:10)
         NERR = NERR + 1
      ENDIF
C
C     TEST CGBSL FOR JOB NOT EQUAL TO 0
C
      CALL CGBSL(AT,LDA,N,ML,MU,IPVT,BT,1)
      INDX = 0
      DO 90 I=1,N
         IF (DELX(C(I),BT(I)) .GT. .0001) INDX=INDX+1
   90 CONTINUE
C
      IF (INDX .NE. 0) THEN
         WRITE (LUN,201) KPROG(11:14),KFAIL(12:19)
         NERR = NERR + 1
      ENDIF
C
C     TEST CGBDI
C
      CALL CGBDI(AT,LDA,N,ML,MU,IPVT,DET)
      INDX = 0
      DO 110 I=1,2
         IF (DELX(DC(I),DET(I)) .GT. .0001) INDX=INDX+1
  110 CONTINUE
C
      IF (INDX .NE. 0) THEN
         WRITE (LUN,201) KPROG(16:19),KFAIL(21:31)
         NERR = NERR + 1
      ENDIF
C
      IF (KPRINT.GE.2 .OR. NERR.NE.0) WRITE (LUN,200) NERR
      RETURN
C
  200 FORMAT(/' * CGBQC - TEST FOR CGBFA, CGBCO, CGBSL AND CGBDI FOUND '
     1   , I1, ' ERRORS.'/)
  201 FORMAT (/' *** C', A, ' FAILURE - ERROR IN ', A)
      END