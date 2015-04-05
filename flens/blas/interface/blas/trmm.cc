#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(strmm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const float     *ALPHA,
            const float     *A_,
            const INTEGER   *LDA,
            float           *B_,
            const INTEGER   *LDB)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: strmm");

#   ifdef TEST_DIRECT_CBLAS

    char    SIDE_   = toupper(*SIDE);
    char    UPLO_   = toupper(*UPLO);
    char    TRANSA_ = toupper(*TRANSA);
    char    DIAG_   = toupper(*DIAG);

    Side           side   = Side(SIDE_);
    StorageUpLo    upLo   = StorageUpLo(UPLO_);
    Transpose      transA = convertTo<Transpose>(TRANSA_);
    Diag           diag   = Diag(DIAG_);

        cblas_strmm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(diag),
                    *M, *N, *ALPHA,
                    A_, *LDA, B_, *LDB);

#   else
        using std::abs;
        using std::max;

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);
        char    TRANSA_ = toupper(*TRANSA);
        char    DIAG_   = toupper(*DIAG);
        INTEGER nRowA   = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (SIDE_!='L' && SIDE_!='R') {
                info = 1;
            } else if (UPLO_!='U' && UPLO_!='L') {
                info = 2;
            } else if (TRANSA_!='N' && TRANSA_!='T' && TRANSA_!='C') {
                info = 3;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 4;
            } else if (*M<0) {
                info = 5;
            } else if (*N<0) {
                info = 6;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDB<max(INTEGER(1),*M)) {
                info = 11;
            }
            if (info!=0) {
                BLAS(xerbla)("STRMM ", &info);
                return;
            }
#       endif

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      transA = convertTo<Transpose>(TRANSA_);
        Diag           diag   = Diag(DIAG_);
        const INTEGER  ka     = side==Left ? *M : *N;

        STrMatrixConstView  A(SFullConstView(ka, ka, *LDA, A_), upLo, diag);
        SGeMatrixView       B = SFullView(*M, *N, *LDB, B_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            if (side==Left) {
                if (transA==NoTrans) {
                    B = alpha*A*B;
                } else {
                    B = alpha*transpose(A)*B;
                }
            } else if (side==Right) {
                if (transA==NoTrans) {
                    B = alpha*B*A;
                } else {
                    B = alpha*B*transpose(A);
                }
            }
#       else
            blas::mm(side, transA, *ALPHA, A, B);
#       endif
#   endif
}

void
BLAS(dtrmm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const double    *ALPHA,
            const double    *A_,
            const INTEGER   *LDA,
            double          *B_,
            const INTEGER   *LDB)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dtrmm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);
        char    TRANSA_ = toupper(*TRANSA);
        char    DIAG_   = toupper(*DIAG);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      transA = convertTo<Transpose>(TRANSA_);
        Diag           diag   = Diag(DIAG_);

        cblas_dtrmm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(diag),
                    *M, *N, *ALPHA,
                    A_, *LDA, B_, *LDB);

#   else

        using std::abs;
        using std::max;

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);
        char    TRANSA_ = toupper(*TRANSA);
        char    DIAG_   = toupper(*DIAG);
        INTEGER nRowA   = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (SIDE_!='L' && SIDE_!='R') {
                info = 1;
            } else if (UPLO_!='U' && UPLO_!='L') {
                info = 2;
            } else if (TRANSA_!='N' && TRANSA_!='T' && TRANSA_!='C') {
                info = 3;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 4;
            } else if (*M<0) {
                info = 5;
            } else if (*N<0) {
                info = 6;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDB<max(INTEGER(1),*M)) {
                info = 11;
            }
            if (info!=0) {
                BLAS(xerbla)("DTRMM ", &info);
                return;
            }
#       endif

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      transA = convertTo<Transpose>(TRANSA_);
        Diag           diag   = Diag(DIAG_);
        const INTEGER  ka     = side==Left ? *M : *N;

        DTrMatrixConstView  A(DFullConstView(ka, ka, *LDA, A_), upLo, diag);
        DGeMatrixView       B = DFullView(*M, *N, *LDB, B_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            if (side==Left) {
                if (transA==NoTrans) {
                    B = alpha*A*B;
                } else {
                    B = alpha*transpose(A)*B;
                }
            } else if (side==Right) {
                if (transA==NoTrans) {
                    B = alpha*B*A;
                } else {
                    B = alpha*B*transpose(A);
                }
            }
#       else
            blas::mm(side, transA, *ALPHA, A, B);
#       endif
#   endif
}

void
BLAS(ctrmm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *A_,
            const INTEGER   *LDA,
            cfloat          *B_,
            const INTEGER   *LDB)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: ctrmm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);
        char    TRANSA_ = toupper(*TRANSA);
        char    DIAG_   = toupper(*DIAG);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      transA = convertTo<Transpose>(TRANSA_);
        Diag           diag   = Diag(DIAG_);

        cblas_ctrmm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(diag),
                    *M, *N,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(A_), *LDA,
                    reinterpret_cast<float *>(B_), *LDB);

#   else

        using std::abs;
        using std::max;

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);
        char    TRANSA_ = toupper(*TRANSA);
        char    DIAG_   = toupper(*DIAG);
        INTEGER nRowA   = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (SIDE_!='L' && SIDE_!='R') {
                info = 1;
            } else if (UPLO_!='U' && UPLO_!='L') {
                info = 2;
            } else if (TRANSA_!='N' && TRANSA_!='T' && TRANSA_!='C') {
                info = 3;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 4;
            } else if (*M<0) {
                info = 5;
            } else if (*N<0) {
                info = 6;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDB<max(INTEGER(1),*M)) {
                info = 11;
            }
            if (info!=0) {
                BLAS(xerbla)("CTRMM ", &info);
                return;
            }
#       endif

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      transA = convertTo<Transpose>(TRANSA_);
        Diag           diag   = Diag(DIAG_);
        const INTEGER  ka     = side==Left ? *M : *N;

        CTrMatrixConstView  A(CFullConstView(ka, ka, *LDA, A_), upLo, diag);
        CGeMatrixView       B = CFullView(*M, *N, *LDB, B_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            if (side==Left) {
                if (transA==NoTrans) {
                    B = alpha*A*B;
                } else if (transA==Trans) {
                    B = alpha*transpose(A)*B;
                } else if (transA==ConjTrans) {
                    B = alpha*conjTrans(A)*B;
                }
            } else if (side==Right) {
                if (transA==NoTrans) {
                    B = alpha*B*A;
                } else if (transA==Trans) {
                    B = alpha*B*transpose(A);
                } else if (transA==ConjTrans) {
                    B = alpha*B*conjTrans(A);
                }
            }
#       else
            blas::mm(side, transA, *ALPHA, A, B);
#       endif
#   endif
}

void
BLAS(ztrmm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *A_,
            const INTEGER   *LDA,
            cdouble         *B_,
            const INTEGER   *LDB)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: ztrmm");

#   ifdef TEST_DIRECT_CBLAS

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);
        char    TRANSA_ = toupper(*TRANSA);
        char    DIAG_   = toupper(*DIAG);

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      transA = convertTo<Transpose>(TRANSA_);
        Diag           diag   = Diag(DIAG_);

        cblas_ztrmm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(side),
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(diag),
                    *M, *N,
                    reinterpret_cast<const double *>(ALPHA),
                    reinterpret_cast<const double *>(A_), *LDA,
                    reinterpret_cast<double *>(B_), *LDB);

#   else

        using std::abs;
        using std::max;

        char    SIDE_   = toupper(*SIDE);
        char    UPLO_   = toupper(*UPLO);
        char    TRANSA_ = toupper(*TRANSA);
        char    DIAG_   = toupper(*DIAG);
        INTEGER nRowA   = (SIDE_=='L') ? *M : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (SIDE_!='L' && SIDE_!='R') {
                info = 1;
            } else if (UPLO_!='U' && UPLO_!='L') {
                info = 2;
            } else if (TRANSA_!='N' && TRANSA_!='T' && TRANSA_!='C') {
                info = 3;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 4;
            } else if (*M<0) {
                info = 5;
            } else if (*N<0) {
                info = 6;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDB<max(INTEGER(1),*M)) {
                info = 11;
            }
            if (info!=0) {
                BLAS(xerbla)("ZTRMM ", &info);
                return;
            }
#       endif

        Side           side   = Side(SIDE_);
        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      transA = convertTo<Transpose>(TRANSA_);
        Diag           diag   = Diag(DIAG_);
        const INTEGER  ka     = side==Left ? *M : *N;

        ZTrMatrixConstView  A(ZFullConstView(ka, ka, *LDA, A_), upLo, diag);
        ZGeMatrixView       B = ZFullView(*M, *N, *LDB, B_);

#   ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;

            if (side==Left) {
                if (transA==NoTrans) {
                    B = alpha*A*B;
                } else if (transA==Trans) {
                    B = alpha*transpose(A)*B;
                } else if (transA==ConjTrans) {
                    B = alpha*conjTrans(A)*B;
                }
            } else if (side==Right) {
                if (transA==NoTrans) {
                    B = alpha*B*A;
                } else if (transA==Trans) {
                    B = alpha*B*transpose(A);
                } else if (transA==ConjTrans) {
                    B = alpha*B*conjTrans(A);
                }
            }
#       else
            blas::mm(side, transA, *ALPHA, A, B);
#       endif
#   endif
}


} // extern "C"
