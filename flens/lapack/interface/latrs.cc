#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LATRS         slatrs_
#   define LATRS_REF     slatrs
#   define LATRS_NAME    "SLATRS"
#elif DOUBLE
#   define LATRS         dlatrs_
#   define LATRS_REF     dlatrs
#   define LATRS_NAME    "DLATRS"
#elif COMPLEX_SINGLE
#   define LATRS         clatrs_
#   define LATRS_REF     clatrs
#   define LATRS_NAME    "CLATRS"
#elif COMPLEX_DOUBLE
#   define LATRS         zlatrs_
#   define LATRS_REF     zlatrs
#   define LATRS_NAME    "ZLATRS"
#endif

using namespace flens;

extern "C" {


void
LATRS_REF(char *UPLO, char *TRANS, char *DIAG, char *NORMIN,
          INT *N, FLOAT *A, INT *LDA, FLOAT *X, FLOAT *SCALE,
          FLOAT *CNORM, INT *INFO);


void
LATRS(char *UPLO, char *TRANS, char *DIAG, char *NORMIN,
      INT *N, FLOAT *A, INT *LDA, FLOAT *X, FLOAT *SCALE,
      FLOAT *CNORM, INT *INFO)
{
#   ifndef USE_LATRS_REF
    ASSERT(0);  // not yet implemented
#   else
    LATRS_REF(UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, CNORM, INFO);
#   endif
}

} // extern "C"

