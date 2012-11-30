      DOUBLE PRECISION FUNCTION DSDOT( N, X, INCX, Y, INCY )

      INTEGER               INCX, INCY, N
      REAL                  X(*), Y(*)

      DOUBLE PRECISION      ZERO
      PARAMETER             ( ZERO = 0.0D0 )

      DOUBLE PRECISION      TEMP

      IF( N.GT.0 ) THEN
         CALL DSDOT_SUB( N, X, INCX, Y, INCY, TEMP )
         DSDOT = TEMP
      ELSE
         DSDOT = ZERO
      END IF

      RETURN
      END
