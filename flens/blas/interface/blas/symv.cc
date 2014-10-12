#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssymv)(const char      *UPLO,
            const INTEGER   *N,
            const float     *ALPHA,
            const float     *A_,
            const INTEGER   *LDA,
            const float     *X,
            const INTEGER   *INCX,
            const float     *BETA,
            float           *Y,
            const INTEGER   *INCY)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_ssymv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N, *ALPHA,
                    A_, *LDA, X, *INCX,
                    *BETA, Y, *INCY);

#   else

        using std::abs;
        using std::max;

        char    UPLO_ = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
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
                BLAS(xerbla)("SSYMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        SSyMatrixConstView    A(SFullConstView(*N, *N, A_, *LDA), upLo);
        SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorView      y(SArrayView(*N, Y, abs(*INCY)), *INCY<0);

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
BLAS(dsymv)(const char      *UPLO,
            const INTEGER   *N,
            const double    *ALPHA,
            const double    *A_,
            const INTEGER   *LDA,
            const double    *X,
            const INTEGER   *INCX,
            const double    *BETA,
            double          *Y,
            const INTEGER   *INCY)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_dsymv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N, *ALPHA,
                    A_, *LDA, X, *INCX,
                    *BETA, Y, *INCY);

#   else
        using std::abs;
        using std::max;

        char    UPLO_ = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
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
                BLAS(xerbla)("DSYMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        DSyMatrixConstView    A(DFullConstView(*N, *N, A_, *LDA), upLo);
        DDenseVectorConstView x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorView      y(DArrayView(*N, Y, abs(*INCY)), *INCY<0);

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
