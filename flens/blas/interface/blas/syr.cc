#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyr)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const float     *X,
           const INTEGER   *INCX,
           float           *_A,
           const INTEGER   *LDA)
{
    using std::abs;
    using std::max;

    INTEGER info  = 0;
    char    _UPLO = toupper(*UPLO);

    if (_UPLO!='U' && _UPLO!='L') {
        info = 1;
    } else if (*N<0) {
        info = 2;
    } else if (*INCX==0) {
        info = 5;
    } else if (*LDA<max(INTEGER(1),*N)) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("SSYR  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    SSyMatrixView         A(SFullView(*N, *N, _A, *LDA), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;

    A += alpha*x*transpose(x);
#   else
    blas::r(*ALPHA, x, A);
#   endif
}

void
BLAS(dsyr)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const double    *X,
           const INTEGER   *INCX,
           double          *_A,
           const INTEGER   *LDA)
{
    using std::abs;
    using std::max;

    INTEGER info  = 0;
    char    _UPLO = toupper(*UPLO);

    if (_UPLO!='U' && _UPLO!='L') {
        info = 1;
    } else if (*N<0) {
        info = 2;
    } else if (*INCX==0) {
        info = 5;
    } else if (*LDA<max(INTEGER(1),*N)) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("DSYR  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    DDenseVectorConstView x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    DSyMatrixView         A(DFullView(*N, *N, _A, *LDA), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;

    A += alpha*x*transpose(x);
#   else
    blas::r(*ALPHA, x, A);
#   endif
}

} // extern "C"
