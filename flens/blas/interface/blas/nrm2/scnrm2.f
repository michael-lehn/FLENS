      REAL FUNCTION SCNRM2( N, X, INCX )

      INTEGER       INCX, N
      COMPLEX       X(*)

      REAL          TEMP

      CALL SCNRM2_SUB( N, X, INCX, TEMP )
      SCNRM2 = TEMP

      RETURN
      END
