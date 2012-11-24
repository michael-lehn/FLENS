      INTEGER FUNCTION IZAMAX( N, X, INCX )

      INTEGER           INCX, N
      DOUBLE COMPLEX    X(*)

      INTEGER           TEMP

      CALL IZAMAX_SUB( N, X, INCX, TEMP )
      IZAMAX = TEMP

      RETURN
      END
