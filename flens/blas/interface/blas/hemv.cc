#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chemv)(const char      *UPLO,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *_A,
            const INTEGER   *LDA,
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
    } else if (*LDA<max(INTEGER(1),*N)) {
        info = 5;
    } else if (*INCX==0) {
        info = 7;
    } else if (*INCY==0) {
        info = 10;
    }
    if (info!=0) {
        BLAS(xerbla)("CHEMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    CHeMatrixConstView    A(CFullConstView(*N, *N, _A, *LDA), upLo);
    CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorView      y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

    // if you only want to test FLENS-BLAS just call
    blas::mv(*ALPHA, A, x, *BETA, y);
}

void
BLAS(zhemv)(const char      *UPLO,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *_A,
            const INTEGER   *LDA,
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
    } else if (*LDA<max(INTEGER(1),*N)) {
        info = 5;
    } else if (*INCX==0) {
        info = 7;
    } else if (*INCY==0) {
        info = 10;
    }
    if (info!=0) {
        BLAS(xerbla)("ZHEMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    ZHeMatrixConstView    A(ZFullConstView(*N, *N, _A, *LDA), upLo);
    ZDenseVectorConstView x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView      y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

    // if you only want to test FLENS-BLAS just call
    blas::mv(*ALPHA, A, x, *BETA, y);
}

} // extern "C"
