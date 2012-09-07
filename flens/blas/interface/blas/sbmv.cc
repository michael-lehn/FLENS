#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const float     *_A,
            const INTEGER   *LDA,
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
    } else if (*K<0) {
        info = 3;
    } else if (*LDA<(*K + 1)) {
        info = 6;
    } else if (*INCX==0) {
        info = 8;
    } else if (*INCY==0) {
        info = 11;
    }
    if (info!=0) {
        BLAS(xerbla)("SSBMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::sbmv(ColMajor, upLo,
                  *N, *K,
                  *ALPHA,
                  _A, *LDA,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
}

void
BLAS(dsbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const double    *_A,
            const INTEGER   *LDA,
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
    } else if (*K<0) {
        info = 3;
    } else if (*LDA<(*K + 1)) {
        info = 6;
    } else if (*INCX==0) {
        info = 8;
    } else if (*INCY==0) {
        info = 11;
    }
    if (info!=0) {
        BLAS(xerbla)("DSBMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::sbmv(ColMajor, upLo,
                  *N, *K,
                  *ALPHA,
                  _A, *LDA,
                  X, *INCX,
                  *BETA,
                  Y, *INCY);
}

} // extern "C"
