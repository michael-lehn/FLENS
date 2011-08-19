#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GEEQU         sgeequ_
#   define GEEQU_REF     sgeequ
#   define GEEQU_NAME    "SGEEQU"
#elif DOUBLE
#   define GEEQU         dgeequ_
#   define GEEQU_REF     dgeequ
#   define GEEQU_NAME    "DGEEQU"
#elif COMPLEX_SINGLE
#   define GEEQU         cgeequ_
#   define GEEQU_REF     cgeequ
#   define GEEQU_NAME    "CGEEQU"
#elif COMPLEX_DOUBLE
#   define GEEQU         zgeequ_
#   define GEEQU_REF     zgeequ
#   define GEEQU_NAME    "ZGEEQU"
#endif

using namespace flens;

extern "C" {

void
GEEQU_REF(INT *M, INT *N, FLOAT *A, INT *LDA,
          FLOAT *R, FLOAT *C, FLOAT *ROWCND, FLOAT *COLCND, FLOAT *AMAX,
          INT *INFO);


void
geequErrorCheck(INT *M, INT *N, FLOAT *A, INT *LDA,
                FLOAT *R, FLOAT *C, FLOAT *ROWCND, FLOAT *COLCND, FLOAT *AMAX,
                INT *INFO)
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
GEEQU(INT *M, INT *N, FLOAT *A, INT *LDA,
      FLOAT *R, FLOAT *C, FLOAT *ROWCND, FLOAT *COLCND, FLOAT *AMAX,
      INT *INFO)
{
    geequErrorCheck(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GEEQU_NAME, INFO, strlen(GEEQU_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GEEQU_REF
    ASSERT(0);  // not yet implemented
#   else
    GEEQU_REF(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO);
#   endif
}

} // extern "C"

