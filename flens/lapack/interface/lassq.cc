#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LASSQ         slassq_
#   define LASSQ_REF     slassq
#   define LASSQ_NAME    "SLASSQ"
#elif DOUBLE
#   define LASSQ         dlassq_
#   define LASSQ_REF     dlassq
#   define LASSQ_NAME    "DLASSQ"
#elif COMPLEX_SINGLE
#   define LASSQ         classq_
#   define LASSQ_REF     classq
#   define LASSQ_NAME    "CLASSQ"
#elif COMPLEX_DOUBLE
#   define LASSQ         zlassq_
#   define LASSQ_REF     zlassq
#   define LASSQ_NAME    "ZLASSQ"
#endif

using namespace flens;

extern "C" {


void
LASSQ_REF(char *UPLO, INT *M, INT *N, FLOAT *ALPHA, FLOAT *BETA,
          FLOAT *A, INT *LDA);


void
LASSQ(char *UPLO, INT *M, INT *N, FLOAT *ALPHA, FLOAT *BETA,
      FLOAT *A, INT *LDA)
{
#   ifndef USE_LASSQ_REF
    ASSERT(0);  // not yet implemented
#   else
    LASSQ_REF(UPLO, M, N, ALPHA, BETA, A, LDA);
#   endif
}

} // extern "C"

