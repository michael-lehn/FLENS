#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chpmv)(const char      *UPLO,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *AP,
            const cfloat    *X,
            const INTEGER   *INCX,
            const cfloat    *BETA,
            cfloat          *Y,
            const INTEGER   *INCY)
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
        info = 6;
    } else if (*INCY==0) {
        info = 9;
    }
    if (info!=0) {
        BLAS(xerbla)("CHPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    CHpMatrixConstView    A(CPackedConstView(*N, AP), upLo);
    CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorView      y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::mv(*ALPHA, A, x, *BETA, y);
}

void
BLAS(zhpmv)(const char      *UPLO,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *AP,
            const cdouble   *X,
            const INTEGER   *INCX,
            const cdouble   *BETA,
            cdouble         *Y,
            const INTEGER   *INCY)
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
        info = 6;
    } else if (*INCY==0) {
        info = 9;
    }
    if (info!=0) {
        BLAS(xerbla)("ZHPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    ZHpMatrixConstView    A(ZPackedConstView(*N, AP), upLo);
    ZDenseVectorConstView x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView      y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::mv(*ALPHA, A, x, *BETA, y);
}

} // extern "C"
