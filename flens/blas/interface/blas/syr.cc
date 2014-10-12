#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyr)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const float     *X,
           const INTEGER   *INCX,
           float           *A_,
           const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_ssyr(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N, *ALPHA,
                   X, *INCX, A_, *LDA);

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
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("SSYR  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        SDenseVectorConstView x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SSyMatrixView         A(SFullView(*N, *N, A_, *LDA), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*transpose(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

void
BLAS(dsyr)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const double    *X,
           const INTEGER   *INCX,
           double          *A_,
           const INTEGER   *LDA)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_dsyr(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N, *ALPHA,
                   X, *INCX, A_, *LDA);

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
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("DSYR  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        DDenseVectorConstView x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        DSyMatrixView         A(DFullView(*N, *N, A_, *LDA), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*transpose(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

} // extern "C"
