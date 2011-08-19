#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LASET         slaset_
#   define LASET_REF     slaset
#   define LASET_NAME    "SLASET"
#elif DOUBLE
#   define LASET         dlaset_
#   define LASET_REF     dlaset
#   define LASET_NAME    "DLASET"
#elif COMPLEX_SINGLE
#   define LASET         claset_
#   define LASET_REF     claset
#   define LASET_NAME    "CLASET"
#elif COMPLEX_DOUBLE
#   define LASET         zlaset_
#   define LASET_REF     zlaset
#   define LASET_NAME    "ZLASET"
#endif

using namespace flens;

extern "C" {


void
LASET_REF(char *UPLO, INT *M, INT *N, FLOAT *ALPHA, FLOAT *BETA,
          FLOAT *A, INT *LDA);


void
LASET(char *UPLO, INT *M, INT *N, FLOAT *ALPHA, FLOAT *BETA,
      FLOAT *A, INT *LDA)
{
#   ifndef USE_LASET_REF
    ASSERT(0);  // not yet implemented
#   else
    LASET_REF(UPLO, M, N, ALPHA, BETA, A, LDA);
#   endif
}

} // extern "C"

