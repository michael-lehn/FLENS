      DOUBLE PRECISION FUNCTION DDOT( N, X, INCX, Y, INCY )

      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   X( * ), Y( * )

      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )

      DOUBLE PRECISION   TEMP
      EXTERNAL           DDOT_SUB

      IF( N.GT.0 ) THEN
         CALL DDOT_SUB( N, X, INCX, Y, INCY, TEMP )
         DDOT = TEMP
      ELSE
         DDOT = ZERO
      END IF

      RETURN
      END
