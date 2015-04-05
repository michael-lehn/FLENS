#include <flens/blas/interface/blas/config.h>

using namespace flens;

extern "C" {

void
BLAS(cher)(const char      *UPLO,
           const INTEGER   *N,
           const float     *ALPHA,
           const cfloat    *X,
           const INTEGER   *INCX,
           cfloat          *A_,
           const INTEGER   *LDA)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: cher");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_cher(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N,
                   *ALPHA,
                   reinterpret_cast<const float *>(X), *INCX,
                   reinterpret_cast<float *>(A_), *LDA);

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
                BLAS(xerbla)("CHER  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        CDenseVectorConstView x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        CHeMatrixView         A(CFullView(*N, *N, *LDA, A_), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

void
BLAS(zher)(const char      *UPLO,
           const INTEGER   *N,
           const double    *ALPHA,
           const cdouble   *X,
           const INTEGER   *INCX,
           cdouble         *A_,
           const INTEGER   *LDA)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zher");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_zher(CBLAS_ORDER::CblasColMajor,
                   cxxblas::CBLAS::getCblasType(upLo),
                   *N,
                   *ALPHA,
                   reinterpret_cast<const double *>(X), *INCX,
                   reinterpret_cast<double *>(A_), *LDA);

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
                BLAS(xerbla)("ZHER  ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);

        ZDenseVectorConstView x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);
        ZHeMatrixView         A(ZFullView(*N, *N, *LDA, A_), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            A += alpha*x*conjTrans(x);
#       else
            blas::r(*ALPHA, x, A);
#       endif
#   endif
}

} // extern "C"
