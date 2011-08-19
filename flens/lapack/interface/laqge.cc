#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LAQGE         slaqge_
#   define LAQGE_REF     slaqge
#   define LAQGE_NAME    "SLAQGE"
#elif DOUBLE
#   define LAQGE         dlaqge_
#   define LAQGE_REF     dlaqge
#   define LAQGE_NAME    "DLAQGE"
#elif COMPLEX_SINGLE
#   define LAQGE         claqge_
#   define LAQGE_REF     claqge
#   define LAQGE_NAME    "CLAQGE"
#elif COMPLEX_DOUBLE
#   define LAQGE         zlaqge_
#   define LAQGE_REF     zlaqge
#   define LAQGE_NAME    "ZLAQGE"
#endif

using namespace flens;

extern "C" {


void
LAQGE_REF(char *NORM, char *UPLO, char *DIAG,
          INT *M, INT *N, FLOAT *A, INT *LDA,
          FLOAT *WORK);


void
LAQGE(char *NORM, char *UPLO, char *DIAG,
      INT *M, INT *N, FLOAT *A, INT *LDA,
      FLOAT *WORK)
{
#   ifndef USE_LAQGE_REF
    ASSERT(0);  // not yet implemented
#   else
    LAQGE_REF(NORM, UPLO, DIAG, M, N, A, LDA, WORK);
#   endif
}

} // extern "C"

