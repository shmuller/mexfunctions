*DECK F4S
      REAL FUNCTION F4S (X)
C***BEGIN PROLOGUE  F4S
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  F4S
      REAL X
C***FIRST EXECUTABLE STATEMENT  F4S
      IF(X.EQ..33E+00) GO TO 10
      F4S = ABS(X-0.33E+00)**(-0.999E+00)
      RETURN
   10 F4S=0.0
      RETURN
      END