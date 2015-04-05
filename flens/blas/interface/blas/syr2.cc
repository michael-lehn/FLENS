#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyr2)(const char      *UPLO,
            const INTEGER   *N,
            const float     *ALPHA,
            const float     *X,
            const INTEGER   *INCX,
            const float     *Y,
            const INTEGER   *INCY,
            float           *A_,
            const INTEGER   *LDA)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: ssyr2");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_ssyr2(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N, *ALPHA,
                   X, *INCX, Y, *INCY, A_, *LDA);

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
                info = 5;
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("SSYR2 ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SDenseVectorConstView y(SConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        SSyMatrixView         A(SFullView(*N, *N, *LDA, A_), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;

            A += alpha*x*transpose(y) + alpha*y*transpose(x);
#       else
            blas::r2(*ALPHA, x, y, A);
#       endif
#   endif
}

void
BLAS(dsyr2)(const char      *UPLO,
            const INTEGER   *N,
            const double    *ALPHA,
            const double    *X,
            const INTEGER   *INCX,
            const double    *Y,
            const INTEGER   *INCY,
            double          *A_,
            const INTEGER   *LDA)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dsyr2");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_dsyr2(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N, *ALPHA,
                   X, *INCX, Y, *INCY, A_, *LDA);

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
                info = 5;
            } else if (*INCY==0) {
                info = 7;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("DSYR2 ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        DDenseVectorConstView x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        DDenseVectorConstView y(DConstArrayView(*N, Y, abs(*INCY)), *INCY<0);
        DSyMatrixView         A(DFullView(*N, *N, *LDA, A_), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;

            A += alpha*x*transpose(y) + alpha*y*transpose(x);
#       else
            blas::r2(*ALPHA, x, y, A);
#       endif
#   endif
}

} // extern "C"
