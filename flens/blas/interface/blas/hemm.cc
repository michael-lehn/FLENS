#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(chemm)(const char      *SIDE,
            const char      *UPLO,
            const INTEGER   *M,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *A_,
            const INTEGER   *LDA,
            const cfloat    *B_,
            const INTEGER   *LDB,
            const cfloat    *BETA,
            cfloat          *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: chemm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_chemm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    *M, *N,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(A_), *LDA,
                    reinterpret_cast<const float *>(B_), *LDB,
                    reinterpret_cast<const float *>(BETA),
                    reinterpret_cast<float *>(C_), *LDC);

#   else

        using std::abs;
        using std::max;

        char    SIDE_ = toupper(*SIDE);
        char    UPLO_ = toupper(*UPLO);
        INTEGER nRowA = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;
            if (SIDE_!='L' && SIDE_!='R') {
                info = 1;
            } else if (UPLO_!='U' && UPLO_!='L') {
                info = 2;
            } else if (*M<0) {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 7;
            } else if (*LDB<max(INTEGER(1),*M)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*M)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("CHEMM ", &info);
                return;
            }
#       endif

        Side           side = Side(SIDE_);
        StorageUpLo    upLo = StorageUpLo(UPLO_);
        const INTEGER  ka   = side==Left ? *M : *N;

        CHeMatrixConstView  A(CFullConstView(ka, ka, *LDA, A_), upLo);
        CGeMatrixConstView  B = CFullConstView(*M, *N, *LDB, B_);
        CGeMatrixView       C = CFullView(*M, *N, *LDC, C_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (side==Left) {
                C = beta*C + alpha*A*B;
            } else if (side==Right) {
                C = beta*C + alpha*B*A;
            }
#       else
            blas::mm(side, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}


void
BLAS(zhemm)(const char      *SIDE,
            const char      *UPLO,
            const INTEGER   *M,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *A_,
            const INTEGER   *LDA,
            const cdouble   *B_,
            const INTEGER   *LDB,
            const cdouble   *BETA,
            cdouble         *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zhemm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_zhemm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    *M, *N,
                    reinterpret_cast<const double *>(ALPHA),
                    reinterpret_cast<const double *>(A_), *LDA,
                    reinterpret_cast<const double *>(B_), *LDB,
                    reinterpret_cast<const double *>(BETA),
                    reinterpret_cast<double *>(C_), *LDC);
#   else

        using std::abs;
        using std::max;

        char    SIDE_ = toupper(*SIDE);
        char    UPLO_ = toupper(*UPLO);
        INTEGER nRowA = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;
            if (SIDE_!='L' && SIDE_!='R') {
                info = 1;
            } else if (UPLO_!='U' && UPLO_!='L') {
                info = 2;
            } else if (*M<0) {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 7;
            } else if (*LDB<max(INTEGER(1),*M)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*M)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("ZHEMM ", &info);
                return;
            }
#       endif

        Side           side = Side(SIDE_);
        StorageUpLo    upLo = StorageUpLo(UPLO_);
        const INTEGER  ka   = side==Left ? *M : *N;

        ZHeMatrixConstView  A(ZFullConstView(ka, ka, *LDA, A_), upLo);
        ZGeMatrixConstView  B = ZFullConstView(*M, *N, *LDB, B_);
        ZGeMatrixView       C = ZFullView(*M, *N, *LDC, C_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (side==Left) {
                C = beta*C + alpha*A*B;
            } else if (side==Right) {
                C = beta*C + alpha*B*A;
            }
#       else
            blas::mm(side, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}

} // extern "C"
