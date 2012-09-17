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
    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    SSbMatrixConstView    A(SBandConstView(*N, *N, KL, KU,  _A, *LDA), upLo);
    SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    SDenseVectorView      y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::mv(*ALPHA, A, x, *BETA, y);
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
    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    DSbMatrixConstView    A(DBandConstView(*N, *N, KL, KU,  _A, *LDA), upLo);
    DDenseVectorConstView x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    DDenseVectorView      y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

    blas::mv(*ALPHA, A, x, *BETA, y);
}

} // extern "C"
