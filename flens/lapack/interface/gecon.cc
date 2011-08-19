#include <flens/lapack/interface/interface.h>
#include <cctype>

#ifdef SINGLE
#   define GECON         sgecon_
#   define GECON_REF     sgecon
#   define GECON_NAME    "SGECON"
#elif DOUBLE
#   define GECON         dgecon_
#   define GECON_REF     dgecon
#   define GECON_NAME    "DGECON"
#elif COMPLEX_SINGLE
#   define GECON         cgecon_
#   define GECON_REF     cgecon
#   define GECON_NAME    "CGECON"
#elif COMPLEX_DOUBLE
#   define GECON         zgecon_
#   define GECON_REF     zgecon
#   define GECON_NAME    "ZGECON"
#endif

using namespace flens;

extern "C" {

void
GECON_REF(char *NORM, INT *N, FLOAT *A, INT *LDA, FLOAT *ANORM,
          FLOAT *RCOND, FLOAT *WORK, INT *IWORK, INT *INFO);


void
geconErrorCheck(char *NORM, INT *N, FLOAT *A, INT *LDA, FLOAT *ANORM,
                FLOAT *RCOND, FLOAT *WORK, INT *IWORK, INT *INFO)
{
    bool oneNorm = (*NORM=='1') || (toupper(*NORM)=='O');
    bool infNorm = (toupper(*NORM)=='I');

    *INFO = 0;

    if ((!oneNorm) && (!infNorm)) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<max(INT(1),*N)) {
        *INFO = -4;
    } else if (*ANORM<0) {
        *INFO = -5;
    }
}


void
GECON(char *NORM, INT *N, FLOAT *A, INT *LDA, FLOAT *ANORM,
      FLOAT *RCOND, FLOAT *WORK, INT *IWORK, INT *INFO)
{
    geconErrorCheck(NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GECON_NAME, INFO, strlen(GECON_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GECON_REF
    ASSERT(0);  // not yet implemented
#   else
    GECON_REF(NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO);
#   endif
}

} // extern "C"

