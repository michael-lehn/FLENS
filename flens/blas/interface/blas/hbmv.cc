#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
            const cfloat    *ALPHA,
            const cfloat    *A_,
            const INTEGER   *LDA,
            const cfloat    *X,
            const INTEGER   *INCX,
            const cfloat    *BETA,
            cfloat          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: chbmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_chbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N, *K,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(A_), *LDA,
                    reinterpret_cast<const float *>(X), *INCX,
                    reinterpret_cast<const float *>(BETA),
                    reinterpret_cast<float *>(Y), *INCY);

#   else

        using std::abs;
        using std::max;

        char    UPLO_ = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;
            if (UPLO_!='U' && UPLO_!='L') {
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
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);
        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        CHbMatrixConstView    A(CBandConstView(*N, *N, KL, KU,  A_, *LDA),
                                upLo);
        CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        CDenseVectorView      y(CArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            y = beta*y + alpha*A*x;
#       else
            blas::mv(*ALPHA, A, x, *BETA, y);
#       endif
#   endif
}

void
BLAS(zhbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
            const cdouble   *ALPHA,
            const cdouble   *A_,
            const INTEGER   *LDA,
            const cdouble   *X,
            const INTEGER   *INCX,
            const cdouble   *BETA,
            cdouble         *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zhbmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_zhbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N, *K,
                    reinterpret_cast<const double *>(ALPHA),
                    reinterpret_cast<const double *>(A_), *LDA,
                    reinterpret_cast<const double *>(X), *INCX,
                    reinterpret_cast<const double *>(BETA),
                    reinterpret_cast<double *>(Y), *INCY);

#   else

        using std::abs;
        using std::max;

        char    UPLO_ = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;
            if (UPLO_!='U' && UPLO_!='L') {
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
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);
        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        ZHbMatrixConstView    A(ZBandConstView(*N, *N, KL, KU,  A_, *LDA),
                                upLo);
        ZDenseVectorConstView x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZDenseVectorView      y(ZArrayView(*N, Y, abs(*INCY)), *INCY<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            y = beta*y + alpha*A*x;
#       else
            blas::mv(*ALPHA, A, x, *BETA, y);
#       endif
#   endif

}

} // extern "C"
