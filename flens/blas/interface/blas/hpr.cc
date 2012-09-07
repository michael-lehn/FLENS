#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chpr)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const cfloat    *X,
           const INTEGER   *INCX,
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
    }
    if (info!=0) {
        BLAS(xerbla)("CHPR  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::hpr(ColMajor, upLo,
                 *N,
                 *ALPHA,
                 X, *INCX,
                 AP);
}

void
BLAS(zhpr)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const cdouble   *X,
           const INTEGER   *INCX,
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
    }
    if (info!=0) {
        BLAS(xerbla)("ZHPR  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::hpr(ColMajor, upLo,
                 *N,
                 *ALPHA,
                 X, *INCX,
                 AP);
}

} // extern "C"
