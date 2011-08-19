#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LANGE         slange_
#   define LANGE_REF     slange
#   define LANGE_NAME    "SLANGE"
#elif DOUBLE
#   define LANGE         dlange_
#   define LANGE_REF     dlange
#   define LANGE_NAME    "DLANGE"
#elif COMPLEX_SINGLE
#   define LANGE         clange_
#   define LANGE_REF     clange
#   define LANGE_NAME    "CLANGE"
#elif COMPLEX_DOUBLE
#   define LANGE         zlange_
#   define LANGE_REF     zlange
#   define LANGE_NAME    "ZLANGE"
#endif

using namespace flens;

extern "C" {


void
LANGE_REF(char *NORM, INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *WORK);

void
LANGE(char *NORM, INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *WORK)
{
#   ifndef USE_LANGE_REF
    ASSERT(0);  // not yet implemented
#   else
    LANGE_REF(NORM, M, N, A, LDA, WORK);
#   endif
}

} // extern "C"

