      DOUBLE COMPLEX FUNCTION ZDOTC( N, X, INCX, Y, INCY )

      INTEGER            INCX, INCY, N
      DOUBLE COMPLEX     X( * ), Y( * )

      DOUBLE COMPLEX     ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )

      DOUBLE COMPLEX     ZTEMP
      EXTERNAL           ZDOTC_SUB

      IF( N.GT.0 ) THEN
         CALL ZDOTC_SUB( N, X, INCX, Y, INCY, ZTEMP )
         ZDOTC = ZTEMP
      ELSE
         ZDOTC = ZERO
      END IF

      RETURN
      END
