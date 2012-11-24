      INTEGER FUNCTION ISAMAX( N, X, INCX )

      INTEGER       INCX, N
      REAL          X(*)

      INTEGER       TEMP

      CALL ISAMAX_SUB( N, X, INCX, TEMP )
      ISAMAX = TEMP

      RETURN
      END
