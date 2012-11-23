      REAL FUNCTION SDOT( N, X, INCX, Y, INCY )

      INTEGER            INCX, INCY, N
      REAL               X( * ), Y( * )

      REAL               ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )

      REAL               STEMP
      EXTERNAL           SDOT_SUB

      IF( N.GT.0 ) THEN
         CALL SDOT_SUB( N, X, INCX, Y, INCY, STEMP )
         SDOT = STEMP
      ELSE
         SDOT = ZERO
      END IF

      RETURN
      END
