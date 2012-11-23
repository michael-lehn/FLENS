      DOUBLE COMPLEX FUNCTION ZDOTU( N, X, INCX, Y, INCY )

      INTEGER            INCX, INCY, N
      DOUBLE COMPLEX     X( * ), Y( * )

      DOUBLE COMPLEX     ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )

      DOUBLE COMPLEX     ZTEMP
      EXTERNAL           ZDOTU_SUB

      IF( N.GT.0 ) THEN
         CALL ZDOTU_SUB( N, X, INCX, Y, INCY, ZTEMP )
         ZDOTU = ZTEMP
      ELSE
         ZDOTU = ZERO
      END IF

      RETURN
      END
