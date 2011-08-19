#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GESVX         sgesvx_
#   define GESVX_REF     sgesvx
#   define GESVX_NAME    "SGESVX"
#elif DOUBLE
#   define GESVX         dgesvx_
#   define GESVX_REF     dgesvx
#   define GESVX_NAME    "DGESVX"
#elif COMPLEX_SINGLE
#   define GESVX         cgesvx_
#   define GESVX_REF     cgesvx
#   define GESVX_NAME    "CGESVX"
#elif COMPLEX_DOUBLE
#   define GESVX         zgesvx_
#   define GESVX_REF     zgesvx
#   define GESVX_NAME    "ZGESVX"
#endif

using namespace flens;

extern "C" {

void
GESVX_REF(char *FACT, char *TRANS, INT *N, INT *NRHS,
          FLOAT *A, INT *LDA, FLOAT *AF, INT *LDAF,
          INT *IPIV,
          char *EQUED,
          FLOAT *R, FLOAT *C,
          FLOAT *B, INT *LDB,
          FLOAT *X, INT *LDX,
          FLOAT *RCOND,
          FLOAT *FERR, FLOAT *BERR,
          FLOAT *WORK, INT *IWORK, INT *INFO);

void
gesvxErrorCheck(char *FACT, char *TRANS, INT *N, INT *NRHS,
                FLOAT *A, INT *LDA, FLOAT *AF, INT *LDAF,
                INT *IPIV,
                char *EQUED,
                FLOAT *R, FLOAT *C,
                FLOAT *B, INT *LDB,
                FLOAT *X, INT *LDX,
                FLOAT *RCOND,
                FLOAT *FERR, FLOAT *BERR,
                FLOAT *WORK, INT *IWORK, INT *INFO)
{
//
//  Test the input parameters.
//
    FLOAT ONE = 1;
    FLOAT SMLNUM, BIGNUM;

    *INFO = 0;
    bool NOFACT = ((*FACT=='N') || (*FACT=='n'));
    bool EQUIL  = ((*FACT=='E') || (*FACT=='e'));
    bool NOTRAN = ((*TRANS=='N') || (*TRANS=='n'));
    bool ROWEQU;
    bool COLEQU;
    if (NOFACT || EQUIL) {
        *EQUED = 'N';
        ROWEQU = false;
        COLEQU = false;
    } else {
        ROWEQU = (*EQUED=='R') || (*EQUED=='r')
              || (*EQUED=='B') || (*EQUED=='b');
        COLEQU = (*EQUED=='C') || (*EQUED=='c')
              || (*EQUED=='B') || (*EQUED=='b');
        SMLNUM = lapack::lamch<FLOAT>(lapack::SafeMin);
        BIGNUM = ONE / SMLNUM;
    }

    if (!NOFACT && !EQUIL && (toupper(*FACT)!='F') ) {
        *INFO = -1;
    } else if (!NOTRAN && (toupper(*TRANS)!='T') && (toupper(*TRANS)!='C')) {
        *INFO = -2;
    } else if (*N<0) {
        *INFO = -3;
    } else if (*NRHS<0) {
        *INFO = -4;
    } else if (*LDA<max(1,*N)) {
        *INFO = -6;
    } else if (*LDAF<max(1,*N)) {
        *INFO = -8;
    } else if ((toupper(*FACT)=='F')
           && !(ROWEQU || COLEQU || toupper(*EQUED)=='N'))
    {
        *INFO = -10;
    } else {
         if (ROWEQU) {
            FLOAT RCMIN = BIGNUM;
            FLOAT RCMAX = FLOAT(0);
            for (INT j=0; j<*N; ++j) {
                RCMIN = min(RCMIN, R[j]);
                RCMAX = max(RCMAX, R[j]);
            }
            if (RCMIN<=FLOAT(0)) {
                *INFO = -11;
            } else if (*N>0) {
            //    ROWCND = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
            } else {
            //    ROWCND = FLOAT(1);
            }
        }
        if (COLEQU && (*INFO==0)) {
            FLOAT RCMIN = BIGNUM;
            FLOAT RCMAX = FLOAT(0);
            for (INT j=0; j<*N; ++j) {
                RCMIN = min(RCMIN, C[j]);
                RCMAX = max(RCMAX, C[j]);
            }
            if (RCMIN<=FLOAT(0)) {
                *INFO = -12;
            } else if (*N>0) {
            //    COLCND = max(RCMIN, SMLNUM) / min(RCMAX, BIGNUM);
            } else {
            //    COLCND = FLOAT(1);
            }
            if (*INFO==0) {
                if (*LDB<max(1, *N)) {
                    *INFO = -14;
                } else if (*LDX<max(1,*N)) {
                    *INFO = -16;
                }
            }
        }
    }

    if (*LDB<max(1, *N)) {
        *INFO = -14;
    }
    if (*LDX<max(1,*N)) {
        *INFO = -16;
    }
}


void
GESVX(char *FACT, char *TRANS, INT *N, INT *NRHS,
      FLOAT *A, INT *LDA, FLOAT *AF, INT *LDAF,
      INT *IPIV,
      char *EQUED,
      FLOAT *R, FLOAT *C,
      FLOAT *B, INT *LDB,
      FLOAT *X, INT *LDX,
      FLOAT *RCOND,
      FLOAT *FERR, FLOAT *BERR,
      FLOAT *WORK, INT *IWORK, INT *INFO)
{
    gesvxErrorCheck(FACT, TRANS, N, NRHS,  A, LDA, AF, LDAF,
                    IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND,
                    FERR, BERR, WORK, IWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GESVX_NAME, INFO, strlen(GESVX_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GESVX_REF
    ASSERT(0);  // not yet implemented
#   else
    GESVX_REF(FACT, TRANS, N, NRHS,  A, LDA, AF, LDAF,
              IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND,
              FERR, BERR, WORK, IWORK, INFO);
#   endif
}

} // extern "C"

