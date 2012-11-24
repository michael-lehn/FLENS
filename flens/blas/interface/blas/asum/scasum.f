      REAL FUNCTION SCASUM( N, X, INCX )

      INTEGER       INCX, N
      COMPLEX       X(*)
      REAL          TEMP


      CALL SCASUM_SUB( N, X, INCX, TEMP )
      SCASUM = TEMP

      RETURN
      END
