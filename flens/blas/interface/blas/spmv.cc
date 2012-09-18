#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sspmv)(const char      *UPLO,
            const INTEGER   *N,
            const float     *ALPHA,
            const float     *AP,
            const float     *X,
            const INTEGER   *INCX,
            const float     *BETA,
            float           *Y,
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
        BLAS(xerbla)("SSPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    SSpMatrixConstView    A(SPackedConstView(*N, AP), upLo);
    SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    SDenseVectorView      y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::mv(*ALPHA, A, x, *BETA, y);
}

void
BLAS(dspmv)(const char      *UPLO,
            const INTEGER   *N,
            const double    *ALPHA,
            const double    *AP,
            const double    *X,
            const INTEGER   *INCX,
            const double    *BETA,
            double          *Y,
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
        BLAS(xerbla)("DSPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    DSpMatrixConstView    A(DPackedConstView(*N, AP), upLo);
    DDenseVectorConstView x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    DDenseVectorView      y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::mv(*ALPHA, A, x, *BETA, y);
}

} // extern "C"
