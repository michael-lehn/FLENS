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

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::spmv(ColMajor, upLo,
                  *N,
                  *ALPHA,
                  AP,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
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

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::spmv(ColMajor, upLo,
                  *N,
                  *ALPHA,
                  AP,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
}

} // extern "C"
