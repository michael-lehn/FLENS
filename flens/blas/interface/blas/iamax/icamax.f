      INTEGER FUNCTION ICAMAX( N, X, INCX )

      INTEGER       INCX, N
      COMPLEX       X(*)

      INTEGER       TEMP

      CALL ICAMAX_SUB( N, X, INCX, TEMP )
      ICAMAX = TEMP

      RETURN
      END
