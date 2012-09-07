#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sspr)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const float     *X,
           const INTEGER   *INCX,
           float           *AP)
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
        BLAS(xerbla)("SSPR  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::spr(ColMajor, upLo,
                 *N,
                 *ALPHA,
                 X, *INCX,
                 AP);
}

void
BLAS(dspr)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const double    *X,
           const INTEGER   *INCX,
           double          *AP)
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
        BLAS(xerbla)("DSPR  ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::spr(ColMajor, upLo,
                 *N,
                 *ALPHA,
                 X, *INCX,
                 AP);
}

} // extern "C"
