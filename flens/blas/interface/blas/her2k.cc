#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(cher2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const cfloat    *ALPHA,
             const cfloat    *A_,
             const INTEGER   *LDA,
             const cfloat    *B_,
             const INTEGER   *LDB,
             const float     *BETA,
             cfloat          *C_,
             const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: cher2k");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_cher2k(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(A_), *LDA,
                    reinterpret_cast<const float *>(B_), *LDB,
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
            } else if (*LDB<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("CHER2K", &info);
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

        CGeMatrixConstView  B = CFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               *LDB,
                                               B_);

        CHeMatrixView       C(CFullView(*N, *N, *LDC, C_), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto alpha_ = conj(alpha);
            const auto beta   = *BETA;

            if (beta==float(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*conjTrans(B) + alpha_*B*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C += alpha*conjTrans(A)*B + alpha_*conjTrans(B)*A;
                }
            } else if (beta==float(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*conjTrans(B) + alpha_*B*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = alpha*conjTrans(A)*B + alpha_*conjTrans(B)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*conjTrans(B) + alpha_*B*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = beta*C + alpha*conjTrans(A)*B + alpha_*conjTrans(B)*A;
                }
            }
#       else
            blas::r2k(trans, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}


void
BLAS(zher2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const cdouble   *ALPHA,
             const cdouble   *A_,
             const INTEGER   *LDA,
             const cdouble   *B_,
             const INTEGER   *LDB,
             const double    *BETA,
             cdouble         *C_,
             const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zher2k");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);

        cblas_zher2k(CBLAS_ORDER::CblasColMajor,
                     cxxblas::CBLAS::getCblasType(upLo),
                     cxxblas::CBLAS::getCblasType(trans),
                     *N, *K,
                     reinterpret_cast<const double *>(ALPHA),
                     reinterpret_cast<const double *>(A_), *LDA,
                     reinterpret_cast<const double *>(B_), *LDB,
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
            } else if (*LDB<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("ZHER2K", &info);
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

        ZGeMatrixConstView  B = ZFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               *LDB,
                                               B_);

        ZHeMatrixView       C(ZFullView(*N, *N, *LDC, C_), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto alpha_ = conj(alpha);
            const auto beta   = *BETA;

            if (beta==double(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*conjTrans(B) + alpha_*B*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C += alpha*conjTrans(A)*B + alpha_*conjTrans(B)*A;
                }
            } else if (beta==double(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*conjTrans(B) + alpha_*B*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = alpha*conjTrans(A)*B + alpha_*conjTrans(B)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*conjTrans(B) + alpha_*B*conjTrans(A);
                } else if (trans==ConjTrans) {
                    C = beta*C + alpha*conjTrans(A)*B + alpha_*conjTrans(B)*A;
                }
            }
#       else
            blas::r2k(trans, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}

} // extern "C"
