#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chpr2)(const char      *UPLO,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *X,
            const INTEGER   *INCX,
            const cfloat    *Y,
            const INTEGER   *INCY,
            cfloat          *AP)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);

    if (_UPLO!='U' && _UPLO!='L') {
        info = 1;
    } else if (*N<0) {
        info = 2;
    } else if (*INCX==0) {
        info = 5;
    } else if (*INCY==0) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("CHPR2 ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorConstView  y(CConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
    CHpMatrixView          A(CPackedView(*N, AP), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto alpha_ = conj(alpha);

    A += alpha*x*conjTrans(y) + alpha_*y*conjTrans(x);
#   else
    blas::r2(*ALPHA, x, y, A);
#   endif
}

void
BLAS(zhpr2)(const char      *UPLO,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *X,
            const INTEGER   *INCX,
            const cdouble   *Y,
            const INTEGER   *INCY,
            cdouble         *AP)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);

    if (_UPLO!='U' && _UPLO!='L') {
        info = 1;
    } else if (*N<0) {
        info = 2;
    } else if (*INCX==0) {
        info = 5;
    } else if (*INCY==0) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("ZHPR2 ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorConstView  y(ZConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
    ZHpMatrixView          A(ZPackedView(*N, AP), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto alpha_ = conj(alpha);

    A += alpha*x*conjTrans(y) + alpha_*y*conjTrans(x);
#   else
    blas::r2(*ALPHA, x, y, A);
#   endif
}

} // extern "C"
