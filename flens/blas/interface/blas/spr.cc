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
    BLAS_DEBUG_OUT("BLAS INTERFACE: sspr");

#   ifdef TEST_DIRECT_CBLAS

    char    UPLO_   = toupper(*UPLO);

    StorageUpLo    upLo   = StorageUpLo(UPLO_);

    cblas_sspr(CBLAS_ORDER::CblasColMajor,
               cxxblas::CBLAS::getCblasType(upLo),
               *N,
               *ALPHA,
               X, *INCX,
               AP);

#   else

        using std::abs;
        using std::max;

        char    UPLO_  = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
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
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        SSpMatrixView          A(SPackedView(*N, AP), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
               const auto alpha = *ALPHA;

               A += alpha*x*transpose(x);
#       else
               blas::r(*ALPHA, x, A);
#       endif
#   endif
}

void
BLAS(dspr)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const double    *X,
           const INTEGER   *INCX,
           double          *AP)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dspr");

#   ifdef TEST_DIRECT_CBLAS

    char    UPLO_   = toupper(*UPLO);

    StorageUpLo    upLo   = StorageUpLo(UPLO_);

    cblas_dspr(CBLAS_ORDER::CblasColMajor,
               cxxblas::CBLAS::getCblasType(upLo),
               *N,
               *ALPHA,
               X, *INCX,
               AP);

#   else

        using std::abs;
        using std::max;

        char    UPLO_  = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
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
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        DSpMatrixView          A(DPackedView(*N, AP), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
               const auto alpha = *ALPHA;

               A += alpha*x*transpose(x);
#       else
               blas::r(*ALPHA, x, A);
#       endif
#   endif
}

} // extern "C"
