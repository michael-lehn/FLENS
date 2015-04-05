#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(cherk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const cfloat    *A_,
            const INTEGER   *LDA,
            const float     *BETA,
            cfloat          *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: cherk");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_cherk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K,
                    *ALPHA,
                    reinterpret_cast<const float *>(A_), *LDA,
                    *BETA,
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
            } else if (TRANS_!='N' && TRANS_!='C') {
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
                BLAS(xerbla)("CHERK ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        CGeMatrixConstView  A = CFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               *LDA,
                                               A_);

        CHeMatrixView       C(CFullView(*N, *N, *LDC, C_), upLo);

    #   ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;


            if (beta==float(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C += alpha*conjTrans(A)*A;
                }
            } else if (beta==float(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = alpha*conjTrans(A)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = beta*C + alpha*conjTrans(A)*A;
                }
            }
#       else
            blas::rk(trans, *ALPHA, A, *BETA, C);
#       endif
#   endif
}

void
BLAS(zherk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const cdouble   *A_,
            const INTEGER   *LDA,
            const double    *BETA,
            cdouble         *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zherk");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_zherk(CBLAS_ORDER::CblasColMajor,
                     cxxblas::CBLAS::getCblasType(upLo),
                     cxxblas::CBLAS::getCblasType(trans),
                     *N, *K,
                     *ALPHA,
                     reinterpret_cast<const double *>(A_), *LDA,
                     *BETA,
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
            } else if (TRANS_!='N' && TRANS_!='C') {
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
                BLAS(xerbla)("ZHERK ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        ZGeMatrixConstView  A = ZFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               *LDA,
                                               A_);

        ZHeMatrixView       C(ZFullView(*N, *N, *LDC, C_), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==double(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C += alpha*conjTrans(A)*A;
                }
            } else if (beta==double(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = alpha*conjTrans(A)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = beta*C + alpha*conjTrans(A)*A;
                }
            }
#       else
            blas::rk(trans, *ALPHA, A, *BETA, C);
#       endif
#   endif
}

} // extern "C"
