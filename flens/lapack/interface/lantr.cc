#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LANTR         slantr_
#   define LANTR_REF     slantr
#   define LANTR_NAME    "SLANTR"
#elif DOUBLE
#   define LANTR         dlantr_
#   define LANTR_REF     dlantr
#   define LANTR_NAME    "DLANTR"
#elif COMPLEX_SINGLE
#   define LANTR         clantr_
#   define LANTR_REF     clantr
#   define LANTR_NAME    "CLANTR"
#elif COMPLEX_DOUBLE
#   define LANTR         zlantr_
#   define LANTR_REF     zlantr
#   define LANTR_NAME    "ZLANTR"
#endif

using namespace flens;

extern "C" {


void
LANTR_REF(char *NORM, char *UPLO, char *DIAG, INT *M, INT *N,
          FLOAT *A, INT *LDA, FLOAT *WORK);

void
LANTR(char *NORM, char *UPLO, char *DIAG, INT *M, INT *N,
      FLOAT *A, INT *LDA, FLOAT *WORK)
{
#   ifndef USE_LANTR_REF
    ASSERT(0);  // not yet implemented
#   else
    LANTR_REF(NORM, UPLO, DIAG, M, N, A, LDA, WORK);
#   endif
}

} // extern "C"

