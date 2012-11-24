      REAL FUNCTION SASUM( N, X, INCX )

      INTEGER       INCX, N
      REAL          X(*)
      REAL          TEMP


      CALL SASUM_SUB( N, X, INCX, TEMP )
      SASUM = TEMP

      RETURN
      END
