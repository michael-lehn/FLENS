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

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::hpr2(ColMajor, upLo,
                  *N,
                  *ALPHA,
                  X, *INCX,
                  Y, *INCY,
                  AP);
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

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::hpr2(ColMajor, upLo,
                  *N,
                  *ALPHA,
                  X, *INCX,
                  Y, *INCY,
                  AP);
}

} // extern "C"
