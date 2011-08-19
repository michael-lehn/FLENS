#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GEQR2P         sgeqr2p_
#   define GEQR2P_REF     sgeqr2p
#   define GEQR2P_NAME    "SGEQR2P"
#elif DOUBLE
#   define GEQR2P         dgeqr2p_
#   define GEQR2P_REF     dgeqr2p
#   define GEQR2P_NAME    "DGEQR2P"
#elif COMPLEX_SINGLE
#   define GEQR2P         cgeqr2p_
#   define GEQR2P_REF     cgeqr2p
#   define GEQR2P_NAME    "CGEQR2P"
#elif COMPLEX_DOUBLE
#   define GEQR2P         zgeqr2p_
#   define GEQR2P_REF     zgeqr2p
#   define GEQR2P_NAME    "ZGEQR2P"
#endif

using namespace flens;

extern "C" {

void
GEQR2P_REF(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
           FLOAT *WORK, INT *LWORK, INT *INFO);


void
geqr2pErrorCheck(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
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
GEQR2P(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
       FLOAT *WORK, INT *LWORK, INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << GEQR2P_NAME << std::endl;
#   endif

    geqr2pErrorCheck(M, N, A, LDA, TAU, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GEQR2P_NAME, INFO, strlen(GEQR2P_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GEQR2P_REF
    ASSERT(0);  // not yet implemented
#   else
    GEQR2P_REF(M, N, A, LDA, TAU, WORK, LWORK, INFO);
#   endif
}

} // extern "C"

