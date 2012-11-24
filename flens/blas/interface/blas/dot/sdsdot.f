      REAL FUNCTION SDSDOT( N, B, X, INCX, Y, INCY )

      REAL          B
      INTEGER       INCX,INCY,N
      REAL          X(*), Y(*)

      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )

      REAL          TEMP

      IF( N.GT.0 ) THEN
         CALL SDSDOT_SUB( N, B, X, INCX, Y, INCY, TEMP )
         SDSDOT = STEMP
      ELSE
         SDSDOT = ZERO
      END IF

      RETURN
      END
