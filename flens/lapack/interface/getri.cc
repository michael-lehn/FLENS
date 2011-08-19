#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GETRI         sgetri_
#   define GETRI_REF     sgetri
#   define GETRI_NAME    "SGETRI"
#elif DOUBLE
#   define GETRI         dgetri_
#   define GETRI_REF     dgetri
#   define GETRI_NAME    "DGETRI"
#elif COMPLEX_SINGLE
#   define GETRI         cgetri_
#   define GETRI_REF     cgetri
#   define GETRI_NAME    "CGETRI"
#elif COMPLEX_DOUBLE
#   define GETRI         zgetri_
#   define GETRI_REF     zgetri
#   define GETRI_NAME    "ZGETRI"
#endif

using namespace flens;

extern "C" {

void
GETRI_REF(INT *N, FLOAT *A, INT *LDA, INT *IPIV,
          FLOAT *WORK, INT *LWORK, INT *INFO);


void
getriErrorCheck(INT *N, FLOAT *A, INT *LDA, INT *IPIV,
                FLOAT *WORK, INT *LWORK, INT *INFO)
{
    bool lQuery = (*LWORK==-1);
    
    if (*N<0) {
        *INFO = -1;
    } else if (*LDA<std::max(INT(1), *N)) {
        *INFO = -3;
    } else if ((*LWORK<std::max(INT(1), *N)) && (!lQuery)) {
        *INFO = -6;
    }
}


void
GETRI(INT *N, FLOAT *A, INT *LDA, INT *IPIV,
      FLOAT *WORK, INT *LWORK, INT *INFO)
{
    getriErrorCheck(N, A, LDA, IPIV, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GETRI_NAME, INFO, strlen(GETRI_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GETRI_REF
    ASSERT(0);  // not yet implemented
#   else
    GETRI_REF(N, A, LDA, IPIV, WORK, LWORK, INFO);
#   endif
}

} // extern "C"

