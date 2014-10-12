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
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_sspmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N,
                    *ALPHA,
                    AP,
                    X, *INCX,
                    *BETA,
                    Y, *INCY);

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
            } else if (*INCX==0) {
                info = 6;
            } else if (*INCY==0) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("SSPMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        SSpMatrixConstView    A(SPackedConstView(*N, AP), upLo);
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
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_dspmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N,
                    *ALPHA,
                    AP,
                    X, *INCX,
                    *BETA,
                    Y, *INCY);

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
            } else if (*INCX==0) {
                info = 6;
            } else if (*INCY==0) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("DSPMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        DSpMatrixConstView    A(DPackedConstView(*N, AP), upLo);
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
