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
#   ifdef TEST_DIRECT_CBLAS

    char    _UPLO   = toupper(*UPLO);

    StorageUpLo    upLo   = StorageUpLo(_UPLO);

    cblas_sspr(CBLAS_ORDER::CblasColMajor,
               cxxblas::CBLAS::getCblasType(upLo),
               *N,
               *ALPHA,
               X, *INCX,
               AP);

#   else

        using std::abs;
        using std::max;

        char    _UPLO  = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
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

        StorageUpLo  upLo = StorageUpLo(_UPLO);

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

#   ifdef TEST_DIRECT_CBLAS

    char    _UPLO   = toupper(*UPLO);

    StorageUpLo    upLo   = StorageUpLo(_UPLO);

    cblas_dspr(CBLAS_ORDER::CblasColMajor,
               cxxblas::CBLAS::getCblasType(upLo),
               *N,
               *ALPHA,
               X, *INCX,
               AP);

#   else

        using std::abs;
        using std::max;

        char    _UPLO  = toupper(*UPLO);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
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

        StorageUpLo  upLo = StorageUpLo(_UPLO);

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
