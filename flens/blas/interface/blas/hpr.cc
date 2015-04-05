#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chpr)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const cfloat    *X,
           const INTEGER   *INCX,
           cfloat          *AP)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: chpr");

#   ifdef TEST_DIRECT_CBLAS

    char    UPLO_   = toupper(*UPLO);

    StorageUpLo    upLo   = StorageUpLo(UPLO_);

    cblas_chpr(CBLAS_ORDER::CblasColMajor,
               cxxblas::CBLAS::getCblasType(upLo),
               *N, *ALPHA,
               reinterpret_cast<const float *>(X), *INCX,
               reinterpret_cast<float *>(AP));

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
                BLAS(xerbla)("CHPR  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        CHpMatrixView          A(CPackedView(*N, AP), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

void
BLAS(zhpr)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const cdouble   *X,
           const INTEGER   *INCX,
           cdouble         *AP)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zhpr");

#   ifdef TEST_DIRECT_CBLAS

    char    UPLO_   = toupper(*UPLO);

    StorageUpLo    upLo   = StorageUpLo(UPLO_);

    cblas_zhpr(CBLAS_ORDER::CblasColMajor,
               cxxblas::CBLAS::getCblasType(upLo),
               *N, *ALPHA,
               reinterpret_cast<const double *>(X), *INCX,
               reinterpret_cast<double *>(AP));

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
                BLAS(xerbla)("ZHPR  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZHpMatrixView          A(ZPackedView(*N, AP), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

} // extern "C"
