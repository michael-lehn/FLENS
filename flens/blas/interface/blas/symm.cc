#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssymm)(const char      *SIDE,
            const char      *UPLO,
            const INTEGER   *M,
            const INTEGER   *N,
            const float     *ALPHA,
            const float     *A_,
            const INTEGER   *LDA,
            const float     *B_,
            const INTEGER   *LDB,
            const float     *BETA,
            float           *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: ssymm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_ssymm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    *M, *N, *ALPHA,
                    A_, *LDA, B_, *LDB,
                    *BETA, C_, *LDC);

#   else

        using std::abs;
        using std::max;

        char    SIDE_ = toupper(*SIDE);
        char    UPLO_ = toupper(*UPLO);
        INTEGER nRowA = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
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
                BLAS(xerbla)("SSYMM ", &info);
                return;
            }
#       endif

        Side           side = Side(SIDE_);
        StorageUpLo    upLo = StorageUpLo(UPLO_);
        const INTEGER  ka   = side==Left ? *M : *N;

        SSyMatrixConstView  A(SFullConstView(ka, ka, *LDA, A_), upLo);
        SGeMatrixConstView  B = SFullConstView(*M, *N, *LDB, B_);
        SGeMatrixView       C = SFullView(*M, *N, *LDC, C_);

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
BLAS(dsymm)(const char      *SIDE,
            const char      *UPLO,
            const INTEGER   *M,
            const INTEGER   *N,
            const double    *ALPHA,
            const double    *A_,
            const INTEGER   *LDA,
            const double    *B_,
            const INTEGER   *LDB,
            const double    *BETA,
            double          *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dsymm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_dsymm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    *M, *N, *ALPHA,
                    A_, *LDA, B_, *LDB,
                    *BETA, C_, *LDC);

#   else

        using std::abs;
        using std::max;

        char    SIDE_ = toupper(*SIDE);
        char    UPLO_ = toupper(*UPLO);
        INTEGER nRowA = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
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
                BLAS(xerbla)("DSYMM ", &info);
                return;
            }
#       endif

        Side           side = Side(SIDE_);
        StorageUpLo    upLo = StorageUpLo(UPLO_);
        const INTEGER  ka   = side==Left ? *M : *N;

        DSyMatrixConstView  A(DFullConstView(ka, ka, *LDA, A_), upLo);
        DGeMatrixConstView  B = DFullConstView(*M, *N, *LDB, B_);
        DGeMatrixView       C = DFullView(*M, *N, *LDC, C_);

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
BLAS(csymm)(const char      *SIDE,
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
    BLAS_DEBUG_OUT("BLAS INTERFACE: csymm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);

        cblas_csymm(CBLAS_ORDER::CblasColMajor,
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
            INTEGER info  = 0;
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
                BLAS(xerbla)("CSYMM ", &info);
                return;
            }
#       endif

        Side           side = Side(SIDE_);
        StorageUpLo    upLo = StorageUpLo(UPLO_);
        const INTEGER  ka   = side==Left ? *M : *N;

        CSyMatrixConstView  A(CFullConstView(ka, ka, *LDA, A_), upLo);
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
BLAS(zsymm)(const char      *SIDE,
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
    BLAS_DEBUG_OUT("BLAS INTERFACE: zsymm");

#   ifdef TEST_DIRECT_CBLAS

    char    SIDE_   = toupper(*SIDE);
    char    UPLO_   = toupper(*UPLO);

    Side           side   = Side(SIDE_);
    StorageUpLo    upLo   = StorageUpLo(UPLO_);

    cblas_zsymm(CBLAS_ORDER::CblasColMajor,
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
            INTEGER info  = 0;
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
                BLAS(xerbla)("ZSYMM ", &info);
                return;
            }
#       endif

        Side           side = Side(SIDE_);
        StorageUpLo    upLo = StorageUpLo(UPLO_);
        const INTEGER  ka   = side==Left ? *M : *N;

        ZSyMatrixConstView  A(ZFullConstView(ka, ka, *LDA, A_), upLo);
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
