*DECK DF1WS
      DOUBLE PRECISION FUNCTION DF1WS (X)
C***BEGIN PROLOGUE  DF1WS
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DF1WS
      DOUBLE PRECISION X
C***FIRST EXECUTABLE STATEMENT  DF1WS
      DF1WS = 0.00D+00
      IF(X-0.33D+00 .NE. 0.00D+00) DF1WS=ABS(X-0.33D+00)**(-0.999D+00)
      RETURN
      END
