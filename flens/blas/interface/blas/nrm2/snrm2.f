      REAL FUNCTION SNRM2( N, X, INCX )

      INTEGER       INCX, N
      REAL          X(*)

      REAL          TEMP

      CALL SNRM2_SUB( N, X, INCX, TEMP )
      SNRM2 = TEMP

      RETURN
      END
