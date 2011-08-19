#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GEQR2         sgeqr2_
#   define GEQR2_REF     sgeqr2
#   define GEQR2_NAME    "SGEQR2"
#elif DOUBLE
#   define GEQR2         dgeqr2_
#   define GEQR2_REF     dgeqr2
#   define GEQR2_NAME    "DGEQR2"
#elif COMPLEX_SINGLE
#   define GEQR2         cgeqr2_
#   define GEQR2_REF     cgeqr2
#   define GEQR2_NAME    "CGEQR2"
#elif COMPLEX_DOUBLE
#   define GEQR2         zgeqr2_
#   define GEQR2_REF     zgeqr2
#   define GEQR2_NAME    "ZGEQR2"
#endif

using namespace flens;

extern "C" {

void
GEQR2_REF(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
          FLOAT *WORK, INT *INFO);


void
geqr2ErrorCheck(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
                FLOAT *WORK, INT *INFO)
{
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -4;
    }
}

void
GEQR2(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
      FLOAT *WORK, INT *LWORK, INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << GEQR2_NAME << std::endl;
#   endif

    geqr2ErrorCheck(M, N, A, LDA, TAU, WORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GEQR2_NAME, INFO, strlen(GEQR2_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GEQR2_REF
    ASSERT(0);  // not yet implemented
#   else
    GEQR2_REF(M, N, A, LDA, TAU, WORK, INFO);
#   endif
}

} // extern "C"

