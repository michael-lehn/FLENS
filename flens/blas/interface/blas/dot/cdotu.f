      COMPLEX FUNCTION CDOTU( N, X, INCX, Y, INCY )

      INTEGER            INCX, INCY, N
      COMPLEX            X( * ), Y( * )

      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )

      COMPLEX            CTEMP
      EXTERNAL           CDOTU_SUB

      IF( N.GT.0 ) THEN
         CALL CDOTU_SUB( N, X, INCX, Y, INCY, CTEMP )
         CDOTU = CTEMP
      ELSE
         CDOTU = ZERO
      END IF

      RETURN
      END
