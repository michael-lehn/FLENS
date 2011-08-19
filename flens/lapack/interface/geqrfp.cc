#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GEQRFP         sgeqrfp_
#   define GEQRFP_REF     sgeqrfp
#   define GEQRFP_NAME    "SGEQRFP"
#elif DOUBLE
#   define GEQRFP         dgeqrfp_
#   define GEQRFP_REF     dgeqrfp
#   define GEQRFP_NAME    "DGEQRFP"
#elif COMPLEX_SINGLE
#   define GEQRFP         cgeqrfp_
#   define GEQRFP_REF     cgeqrfp
#   define GEQRFP_NAME    "CGEQRFP"
#elif COMPLEX_DOUBLE
#   define GEQRFP         zgeqrfp_
#   define GEQRFP_REF     zgeqrfp
#   define GEQRFP_NAME    "ZGEQRFP"
#endif

using namespace flens;

extern "C" {

void
GEQRFP_REF(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
           FLOAT *WORK, INT *LWORK, INT *INFO);


void
geqrfpErrorCheck(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
                 FLOAT *WORK, INT *LWORK, INT *INFO)
{
    bool lQuery = (*LWORK==-1);
    
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -4;
    } else if ((*LWORK<std::max(INT(1), *N)) && (!lQuery)) {
        *INFO = -7;
    }
}


void
GEQRFP(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
       FLOAT *WORK, INT *LWORK, INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << GEQRFP_NAME << std::endl;
#   endif

    geqrfpErrorCheck(M, N, A, LDA, TAU, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GEQRFP_NAME, INFO, strlen(GEQRFP_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GEQRFP_REF
    ASSERT(0);  // not yet implemented
#   else
    GEQRFP_REF(M, N, A, LDA, TAU, WORK, LWORK, INFO);
#   endif
}

} // extern "C"

