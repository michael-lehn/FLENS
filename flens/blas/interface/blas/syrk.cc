#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const float     *A_,
            const INTEGER   *LDA,
            const float     *BETA,
            float           *C_,
            const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_ssyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K, *ALPHA,
                    A_, *LDA, *BETA, C_, *LDC);

#   else

        using std::abs;
        using std::max;

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        INTEGER nRowA  = (TRANS_=='N') ? *N : *K;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
                info = 2;
            } else if (*N<0) {
                info = 3;
            } else if (*K<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 7;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 10;
            }
            if (info!=0) {
                BLAS(xerbla)("SSYRK ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        SGeMatrixConstView  A = SFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               A_,
                                               *LDA);

        SSyMatrixView       C(SFullView(*N, *N, C_, *LDC), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==float(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C += alpha*transpose(A)*A;
                }
            } else if (beta==float(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = alpha*transpose(A)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = beta*C + alpha*transpose(A)*A;
                }
            }
#       else
            blas::rk(trans, *ALPHA, A, *BETA, C);
#       endif
#   endif
}

void
BLAS(dsyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const double    *A_,
            const INTEGER   *LDA,
            const double    *BETA,
            double          *C_,
            const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_dsyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K, *ALPHA,
                    A_, *LDA, *BETA, C_, *LDC);
#   else
        using std::abs;
        using std::max;

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        INTEGER nRowA  = (TRANS_=='N') ? *N : *K;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
                info = 2;
            } else if (*N<0) {
                info = 3;
            } else if (*K<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 7;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 10;
            }
            if (info!=0) {
                BLAS(xerbla)("DSYRK ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        DGeMatrixConstView  A = DFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               A_,
                                               *LDA);

        DSyMatrixView       C(DFullView(*N, *N, C_, *LDC), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==double(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C += alpha*transpose(A)*A;
                }
            } else if (beta==double(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = alpha*transpose(A)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = beta*C + alpha*transpose(A)*A;
                }
            }
#       else
            blas::rk(trans, *ALPHA, A, *BETA, C);
#       endif
#   endif
}

void
BLAS(csyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const cfloat    *ALPHA,
            const cfloat    *A_,
            const INTEGER   *LDA,
            const cfloat    *BETA,
            cfloat          *C_,
            const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_csyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(A_), *LDA,
                    reinterpret_cast<const float *>(BETA),
                    reinterpret_cast<float *>(C_), *LDC);
#   else

        using std::abs;
        using std::max;

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        INTEGER nRowA  = (TRANS_=='N') ? *N : *K;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T') {
                info = 2;
            } else if (*N<0) {
                info = 3;
            } else if (*K<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 7;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 10;
            }
            if (info!=0) {
                BLAS(xerbla)("CSYRK ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        CGeMatrixConstView  A = CFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               A_,
                                               *LDA);

        CSyMatrixView       C(CFullView(*N, *N, C_, *LDC), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==cfloat(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C += alpha*transpose(A)*A;
                }
            } else if (beta==cfloat(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = alpha*transpose(A)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = beta*C + alpha*transpose(A)*A;
                }
            }
#       else
            blas::rk(trans, *ALPHA, A, *BETA, C);
#       endif
#   endif
}

void
BLAS(zsyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const cdouble   *ALPHA,
            const cdouble   *A_,
            const INTEGER   *LDA,
            const cdouble   *BETA,
            cdouble         *C_,
            const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_zsyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K,
                    reinterpret_cast<const double *>(ALPHA),
                    reinterpret_cast<const double *>(A_), *LDA,
                    reinterpret_cast<const double *>(BETA),
                    reinterpret_cast<double *>(C_), *LDC);
#   else
        using std::abs;
        using std::max;

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        INTEGER nRowA  = (TRANS_=='N') ? *N : *K;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;

            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T') {
                info = 2;
            } else if (*N<0) {
                info = 3;
            } else if (*K<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 7;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 10;
            }
            if (info!=0) {
                BLAS(xerbla)("ZSYRK ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        ZGeMatrixConstView  A = ZFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               A_,
                                               *LDA);

        ZSyMatrixView       C(ZFullView(*N, *N, C_, *LDC), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==cdouble(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C += alpha*transpose(A)*A;
                }
            } else if (beta==cdouble(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = alpha*transpose(A)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(A);
                } else if (trans==Trans || trans==ConjTrans) {
                    C = beta*C + alpha*transpose(A)*A;
                }
            }
#       else
            blas::rk(trans, *ALPHA, A, *BETA, C);
#       endif
#   endif
}


} // extern "C"
