      COMPLEX FUNCTION CDOTC( N, X, INCX, Y, INCY )

      INTEGER            INCX, INCY, N
      COMPLEX            X( * ), Y( * )

      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )

      COMPLEX            CTEMP
      EXTERNAL           CDOTC_SUB

      IF( N.GT.0 ) THEN
         CALL CDOTC_SUB( N, X, INCX, Y, INCY, CTEMP )
         CDOTC = CTEMP
      ELSE
         CDOTC = ZERO
      END IF

      RETURN
      END
