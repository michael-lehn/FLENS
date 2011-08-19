#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LASWP         slaswp_
#   define LASWP_REF     slaswp
#   define LASWP_NAME    "SLASWP"
#elif DOUBLE
#   define LASWP         dlaswp_
#   define LASWP_REF     dlaswp
#   define LASWP_NAME    "DLASWP"
#elif COMPLEX_SINGLE
#   define LASWP         claswp_
#   define LASWP_REF     claswp
#   define LASWP_NAME    "CLASWP"
#elif COMPLEX_DOUBLE
#   define LASWP         zlaswp_
#   define LASWP_REF     zlaswp
#   define LASWP_NAME    "ZLASWP"
#endif

using namespace flens;

extern "C" {


void
LASWP_REF(char *UPLO, char *TRANS, char *DIAG, char *NORMIN,
          INT *N, FLOAT *A, INT *LDA, FLOAT *X, FLOAT *SCALE,
          FLOAT *CNORM, INT *INFO);


void
LASWP(char *UPLO, char *TRANS, char *DIAG, char *NORMIN,
      INT *N, FLOAT *A, INT *LDA, FLOAT *X, FLOAT *SCALE,
      FLOAT *CNORM, INT *INFO)
{
#   ifndef USE_LASWP_REF
    ASSERT(0);  // not yet implemented
#   else
    LASWP_REF(UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, CNORM, INFO);
#   endif
}

} // extern "C"

