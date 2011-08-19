#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LACPY         slacpy_
#   define LACPY_REF     slacpy
#   define LACPY_NAME    "SLACPY"
#elif DOUBLE
#   define LACPY         dlacpy_
#   define LACPY_REF     dlacpy
#   define LACPY_NAME    "DLACPY"
#elif COMPLEX_SINGLE
#   define LACPY         clacpy_
#   define LACPY_REF     clacpy
#   define LACPY_NAME    "CLACPY"
#elif COMPLEX_DOUBLE
#   define LACPY         zlacpy_
#   define LACPY_REF     zlacpy
#   define LACPY_NAME    "ZLACPY"
#endif

using namespace flens;

extern "C" {


void
LACPY_REF(char *UPLO, INT *M, INT *N,
          FLOAT *A, INT *LDA,
          FLOAT *B, INT *LDB);


void
LACPY(char *UPLO, INT *M, INT *N,
      FLOAT *A, INT *LDA,
      FLOAT *B, INT *LDB)
{
#   ifndef USE_LACPY_REF
    ASSERT(0);  // not yet implemented
#   else
    LACPY_REF(UPLO, M, N, A, LDA, B, LDB);
#   endif
}

} // extern "C"

