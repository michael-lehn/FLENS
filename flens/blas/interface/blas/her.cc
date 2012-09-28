#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(cher)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const cfloat    *X,
           const INTEGER   *INCX,
           cfloat          *_A,
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
        BLAS(xerbla)("CHER  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CHeMatrixView         A(CFullView(*N, *N, _A, *LDA), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;

    A += alpha*x*conjTrans(x);
#   else
    blas::r(*ALPHA, x, A);
#   endif
}

void
BLAS(zher)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const cdouble   *X,
           const INTEGER   *INCX,
           cdouble         *_A,
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
        BLAS(xerbla)("ZHER  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    ZDenseVectorConstView x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZHeMatrixView         A(ZFullView(*N, *N, _A, *LDA), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;

    A += alpha*x*conjTrans(x);
#   else
    blas::r(*ALPHA, x, A);
#   endif
}

} // extern "C"
