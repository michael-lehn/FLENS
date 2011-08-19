#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define LANSY         slansy_
#   define LANSY_REF     slansy
#   define LANSY_NAME    "SLANSY"
#elif DOUBLE
#   define LANSY         dlansy_
#   define LANSY_REF     dlansy
#   define LANSY_NAME    "DLANSY"
#elif COMPLEX_SINGLE
#   define LANSY         clansy_
#   define LANSY_REF     clansy
#   define LANSY_NAME    "CLANSY"
#elif COMPLEX_DOUBLE
#   define LANSY         zlansy_
#   define LANSY_REF     zlansy
#   define LANSY_NAME    "ZLANSY"
#endif

using namespace flens;

extern "C" {


void
LANSY_REF(char *NORM, char *UPLO, INT *N, FLOAT *A, INT *LDA, FLOAT *WORK);

void
LANSY(char *NORM, char *UPLO, INT *N, FLOAT *A, INT *LDA, FLOAT *WORK)
{
#   ifndef USE_LANSY_REF
    ASSERT(0);  // not yet implemented
#   else
    LANSY_REF(NORM, UPLO, N, A, LDA, WORK);
#   endif
}

} // extern "C"

