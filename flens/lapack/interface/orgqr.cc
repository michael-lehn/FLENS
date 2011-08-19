#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define ORGQR         sorgqr_
#   define ORGQR_REF     sorgqr
#   define ORGQR_NAME    "SORGQR"
#elif DOUBLE
#   define ORGQR         dorgqr_
#   define ORGQR_REF     dorgqr
#   define ORGQR_NAME    "DORGQR"
#elif COMPLEX_SINGLE
#   define ORGQR         corgqr_
#   define ORGQR_REF     corgqr
#   define ORGQR_NAME    "CORGQR"
#elif COMPLEX_DOUBLE
#   define ORGQR         zorgqr_
#   define ORGQR_REF     zorgqr
#   define ORGQR_NAME    "ZORGQR"
#endif

using namespace flens;

extern "C" {

void
ORGQR_REF(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
          FLOAT *WORK, INT *LWORK, INT *INFO);

void
orgqrErrorCheck(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
                FLOAT *WORK, INT *LWORK, INT *INFO)
{
    bool lQuery = (*LWORK==-1);
    
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*K<0) {
        *INFO = -3;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -5;
    } else if ((*LWORK<std::max(INT(1), *N)) && (!lQuery)) {
        *INFO = -8;
    }
}


void
ORGQR(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
      FLOAT *WORK, INT *LWORK, INT *INFO)
{
    orgqrErrorCheck(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(ORGQR_NAME, INFO, strlen(ORGQR_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_ORGQR_REF
    ASSERT(0);  // not yet implemented
#   else
    ORGQR_REF(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);
#   endif
}

} // extern "C"

