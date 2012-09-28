#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
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
        BLAS(xerbla)("CHBMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    CHbMatrixConstView    A(CBandConstView(*N, *N, KL, KU,  _A, *LDA), upLo);
    CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    CDenseVectorView      y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;
    const auto beta  = *BETA;

    y = beta*y + alpha*A*x;
#   else
    blas::mv(*ALPHA, A, x, *BETA, y);
#   endif
}

void
BLAS(zhbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
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
        BLAS(xerbla)("ZHBMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    ZHbMatrixConstView    A(ZBandConstView(*N, *N, KL, KU,  _A, *LDA), upLo);
    ZDenseVectorConstView x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
    ZDenseVectorView      y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;
    const auto beta  = *BETA;

    y = beta*y + alpha*A*x;
#   else
    blas::mv(*ALPHA, A, x, *BETA, y);
#   endif
}

} // extern "C"
