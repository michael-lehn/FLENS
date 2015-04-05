#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const float     *A_,
            const INTEGER   *LDA,
            const float     *X,
            const INTEGER   *INCX,
            const float     *BETA,
            float           *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: ssbmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_ssbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N, *K,
                    *ALPHA,
                    A_, *LDA,
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
                BLAS(xerbla)("SSBMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);
        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        SSbMatrixConstView    A(SBandConstView(*N, *N, KL, KU,  A_, *LDA),
                                upLo);
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
BLAS(dsbmv)(const char      *UPLO,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const double    *A_,
            const INTEGER   *LDA,
            const double    *X,
            const INTEGER   *INCX,
            const double    *BETA,
            double          *Y,
            const INTEGER   *INCY)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dsbmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_dsbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    *N, *K,
                    *ALPHA,
                    A_, *LDA,
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
                BLAS(xerbla)("DSBMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);
        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        DSbMatrixConstView    A(DBandConstView(*N, *N, KL, KU,  A_, *LDA),
                                upLo);
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
