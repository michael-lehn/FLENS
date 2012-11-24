      INTEGER FUNCTION IDAMAX( N, X, INCX )

      INTEGER               INCX, N
      DOUBLE PRECISION      X(*)

      INTEGER               TEMP

      CALL IDAMAX_SUB( N, X, INCX, TEMP )
      IDAMAX = TEMP

      RETURN
      END
