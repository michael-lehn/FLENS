#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const float     *A_,
            const INTEGER   *LDA,
            const float     *B_,
            const INTEGER   *LDB,
            const float     *BETA,
            float           *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: sgemm");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANSA_ = toupper(*TRANSA);
        const char         TRANSB_ = toupper(*TRANSB);
        const Transpose    transA  = convertTo<Transpose>(TRANSA_);
        const Transpose    transB  = convertTo<Transpose>(TRANSB_);

        cblas_sgemm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(transB),
                    *M, *N, *K, *ALPHA,
                    A_, *LDA, B_, *LDB, *BETA, C_, *LDC);

#   else

        using std::abs;
        using std::max;

        char    TRANSA_ = toupper(*TRANSA);
        char    TRANSB_ = toupper(*TRANSB);
        INTEGER nRowA   = (TRANSA_=='N') ? *M : *K;
        INTEGER nRowB   = (TRANSB_=='N') ? *K : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info = 0;
            if (TRANSA_!='N' && TRANSA_!='C' && TRANSA_!='T') {
                info = 1;
            } else if (TRANSB_!='N' && TRANSB_!='C' && TRANSB_!='T') {
                info = 2;
            } else if (*M<0) {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 8;
            } else if (*LDB<max(INTEGER(1),nRowB)) {
                info = 10;
            } else if (*LDC<max(INTEGER(1),*M)) {
                info = 13;
            }
            if (info!=0) {
                BLAS(xerbla)("SGEMM ", &info);
                return;
            }
#       endif

        Transpose  transA = convertTo<Transpose>(TRANSA_);
        Transpose  transB = convertTo<Transpose>(TRANSB_);

        const bool noTransA = (transA==NoTrans || transA==Conj);
        const bool noTransB = (transB==NoTrans || transB==Conj);

        const SGeMatrixConstView  A = SFullConstView(noTransA ? *M : *K,
                                                     noTransA ? *K : *M,
                                                     *LDA,
                                                     A_);

        const SGeMatrixConstView  B = SFullConstView(noTransB ? *K : *N,
                                                     noTransB ? *N : *K,
                                                     *LDB,
                                                     B_);

        SGeMatrixView       C = SFullView(*M, *N, *LDC, C_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (transA==NoTrans && transB==NoTrans) {
                C = beta*C + alpha*A*B;
            } else if (transA==Trans && transB==NoTrans) {
                C = beta*C + alpha*transpose(A)*B;
            } else if (transA==ConjTrans && transB==NoTrans) {
                C = beta*C + alpha*conjTrans(A)*B;
            } else if (transA==NoTrans && transB==Trans) {
                C = beta*C + alpha*A*transpose(B);
            } else if (transA==Trans && transB==Trans) {
                C = beta*C + alpha*transpose(A)*transpose(B);
            } else if (transA==ConjTrans && transB==Trans) {
                C = beta*C + alpha*conjTrans(A)*transpose(B);
            } else if (transA==NoTrans && transB==ConjTrans) {
                C = beta*C + alpha*A*conjTrans(B);
            } else if (transA==Trans && transB==ConjTrans) {
                C = beta*C + alpha*transpose(A)*conjTrans(B);
            } else if (transA==ConjTrans && transB==ConjTrans) {
                C = beta*C + alpha*conjTrans(A)*conjTrans(B);
            }
#       else
            blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}


void
BLAS(dgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const double    *A_,
            const INTEGER   *LDA,
            const double    *B_,
            const INTEGER   *LDB,
            const double    *BETA,
            double          *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dgemm");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANSA_ = toupper(*TRANSA);
        const char         TRANSB_ = toupper(*TRANSB);
        const Transpose    transA  = convertTo<Transpose>(TRANSA_);
        const Transpose    transB  = convertTo<Transpose>(TRANSB_);

        cblas_dgemm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(transB), *M, *N, *K, *ALPHA,
                    A_, *LDA, B_, *LDB, *BETA, C_, *LDC);

#   else

        using std::abs;
        using std::max;


        char    TRANSA_ = toupper(*TRANSA);
        char    TRANSB_ = toupper(*TRANSB);
        INTEGER nRowA   = (TRANSA_=='N') ? *M : *K;
        INTEGER nRowB   = (TRANSB_=='N') ? *K : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;
            if (TRANSA_!='N' && TRANSA_!='C' && TRANSA_!='T') {
                info = 1;
            } else if (TRANSB_!='N' && TRANSB_!='C' && TRANSB_!='T') {
                info = 2;
            } else if (*M<0) {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 8;
            } else if (*LDB<max(INTEGER(1),nRowB)) {
                info = 10;
            } else if (*LDC<max(INTEGER(1),*M)) {
                info = 13;
            }
            if (info!=0) {
                BLAS(xerbla)("DGEMM ", &info);
                return;
            }
#       endif


        Transpose  transA = convertTo<Transpose>(TRANSA_);
        Transpose  transB = convertTo<Transpose>(TRANSB_);

        const bool noTransA = (transA==NoTrans || transA==Conj);
        const bool noTransB = (transB==NoTrans || transB==Conj);

        const DGeMatrixConstView  A = DFullConstView(noTransA ? *M : *K,
                                                     noTransA ? *K : *M,
                                                     *LDA,
                                                     A_);

        const DGeMatrixConstView  B = DFullConstView(noTransB ? *K : *N,
                                                     noTransB ? *N : *K,
                                                     *LDB,
                                                     B_);

        DGeMatrixView  C = DFullView(*M, *N, *LDC, C_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (transA==NoTrans && transB==NoTrans) {
                C = beta*C + alpha*A*B;
            } else if (transA==Trans && transB==NoTrans) {
                C = beta*C + alpha*transpose(A)*B;
            } else if (transA==ConjTrans && transB==NoTrans) {
                C = beta*C + alpha*conjTrans(A)*B;
            } else if (transA==NoTrans && transB==Trans) {
                C = beta*C + alpha*A*transpose(B);
            } else if (transA==Trans && transB==Trans) {
                C = beta*C + alpha*transpose(A)*transpose(B);
            } else if (transA==ConjTrans && transB==Trans) {
                C = beta*C + alpha*conjTrans(A)*transpose(B);
            } else if (transA==NoTrans && transB==ConjTrans) {
                C = beta*C + alpha*A*conjTrans(B);
            } else if (transA==Trans && transB==ConjTrans) {
                C = beta*C + alpha*transpose(A)*conjTrans(B);
            } else if (transA==ConjTrans && transB==ConjTrans) {
                C = beta*C + alpha*conjTrans(A)*conjTrans(B);
            }
#       else
            blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}

void
BLAS(cgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const cfloat    *ALPHA,
            const cfloat    *A_,
            const INTEGER   *LDA,
            const cfloat    *B_,
            const INTEGER   *LDB,
            const cfloat    *BETA,
            cfloat          *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: cgemm");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANSA_ = toupper(*TRANSA);
        const char         TRANSB_ = toupper(*TRANSB);
        const Transpose    transA  = convertTo<Transpose>(TRANSA_);
        const Transpose    transB  = convertTo<Transpose>(TRANSB_);

        cblas_cgemm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(transB),
                    *M, *N, *K,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(A_), *LDA,
                    reinterpret_cast<const float *>(B_), *LDB,
                    reinterpret_cast<const float *>(BETA),
                    reinterpret_cast<float *>(C_), *LDC);

#   else

        using std::abs;
        using std::max;

        char    TRANSA_ = toupper(*TRANSA);
        char    TRANSB_ = toupper(*TRANSB);
        INTEGER nRowA   = (TRANSA_=='N') ? *M : *K;
        INTEGER nRowB   = (TRANSB_=='N') ? *K : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;
            if (TRANSA_!='N' && TRANSA_!='C' && TRANSA_!='T') {
                info = 1;
            } else if (TRANSB_!='N' && TRANSB_!='C' && TRANSB_!='T') {
                info = 2;
            } else if (*M<0) {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 8;
            } else if (*LDB<max(INTEGER(1),nRowB)) {
                info = 10;
            } else if (*LDC<max(INTEGER(1),*M)) {
                info = 13;
            }
            if (info!=0) {
                BLAS(xerbla)("CGEMM ", &info);
                return;
            }
#       endif

        Transpose  transA = convertTo<Transpose>(TRANSA_);
        Transpose  transB = convertTo<Transpose>(TRANSB_);

        const bool noTransA = (transA==NoTrans || transA==Conj);
        const bool noTransB = (transB==NoTrans || transB==Conj);

        const CGeMatrixConstView  A = CFullConstView(noTransA ? *M : *K,
                                                     noTransA ? *K : *M,
                                                     *LDA,
                                                     A_);

        const CGeMatrixConstView  B = CFullConstView(noTransB ? *K : *N,
                                                     noTransB ? *N : *K,
                                                     *LDB,
                                                     B_);

        CGeMatrixView       C = CFullView(*M, *N, *LDC, C_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (transA==NoTrans && transB==NoTrans) {
                C = beta*C + alpha*A*B;
            } else if (transA==Trans && transB==NoTrans) {
                C = beta*C + alpha*transpose(A)*B;
            } else if (transA==ConjTrans && transB==NoTrans) {
                C = beta*C + alpha*conjTrans(A)*B;
            } else if (transA==NoTrans && transB==Trans) {
                C = beta*C + alpha*A*transpose(B);
            } else if (transA==Trans && transB==Trans) {
                C = beta*C + alpha*transpose(A)*transpose(B);
            } else if (transA==ConjTrans && transB==Trans) {
                C = beta*C + alpha*conjTrans(A)*transpose(B);
            } else if (transA==NoTrans && transB==ConjTrans) {
                C = beta*C + alpha*A*conjTrans(B);
            } else if (transA==Trans && transB==ConjTrans) {
                C = beta*C + alpha*transpose(A)*conjTrans(B);
            } else if (transA==ConjTrans && transB==ConjTrans) {
                C = beta*C + alpha*conjTrans(A)*conjTrans(B);
            }
#       else
            blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}


void
BLAS(zgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const cdouble   *ALPHA,
            const cdouble   *A_,
            const INTEGER   *LDA,
            const cdouble   *B_,
            const INTEGER   *LDB,
            const cdouble   *BETA,
            cdouble         *C_,
            const INTEGER   *LDC)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: zgbmv");

#   ifdef TEST_DIRECT_CBLAS

        const char         TRANSA_ = toupper(*TRANSA);
        const char         TRANSB_ = toupper(*TRANSB);
        const Transpose    transA  = convertTo<Transpose>(TRANSA_);
        const Transpose    transB  = convertTo<Transpose>(TRANSB_);

        cblas_zgemm(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(transA),
                    cxxblas::CBLAS::getCblasType(transB),
                    *M, *N, *K,
                    reinterpret_cast<const double *>(ALPHA),
                    reinterpret_cast<const double *>(A_), *LDA,
                    reinterpret_cast<const double *>(B_), *LDB,
                    reinterpret_cast<const double *>(BETA),
                    reinterpret_cast<double *>(C_), *LDC);

#   else

        using std::abs;
        using std::max;


        char    TRANSA_ = toupper(*TRANSA);
        char    TRANSB_ = toupper(*TRANSB);
        INTEGER nRowA   = (TRANSA_=='N') ? *M : *K;
        INTEGER nRowB   = (TRANSB_=='N') ? *K : *N;

#       ifndef NO_INPUT_CHECK
            INTEGER info    = 0;
            if (TRANSA_!='N' && TRANSA_!='C' && TRANSA_!='T') {
                info = 1;
            } else if (TRANSB_!='N' && TRANSB_!='C' && TRANSB_!='T') {
                info = 2;
            } else if (*M<0) {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<max(INTEGER(1),nRowA)) {
                info = 8;
            } else if (*LDB<max(INTEGER(1),nRowB)) {
                info = 10;
            } else if (*LDC<max(INTEGER(1),*M)) {
                info = 13;
            }
            if (info!=0) {
                BLAS(xerbla)("ZGEMM ", &info);
                return;
            }
#       endif



        Transpose  transA = convertTo<Transpose>(TRANSA_);
        Transpose  transB = convertTo<Transpose>(TRANSB_);

        const bool noTransA = (transA==NoTrans || transA==Conj);
        const bool noTransB = (transB==NoTrans || transB==Conj);

        const ZGeMatrixConstView  A = ZFullConstView(noTransA ? *M : *K,
                                                     noTransA ? *K : *M,
                                                     *LDA,
                                                     A_);

        const ZGeMatrixConstView  B = ZFullConstView(noTransB ? *K : *N,
                                                     noTransB ? *N : *K,
                                                     *LDB,
                                                     B_);

        ZGeMatrixView  C = ZFullView(*M, *N, *LDC, C_);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha = *ALPHA;
            const auto beta  = *BETA;

            if (transA==NoTrans && transB==NoTrans) {
                C = beta*C + alpha*A*B;
            } else if (transA==Trans && transB==NoTrans) {
                C = beta*C + alpha*transpose(A)*B;
            } else if (transA==ConjTrans && transB==NoTrans) {
                C = beta*C + alpha*conjTrans(A)*B;
            } else if (transA==NoTrans && transB==Trans) {
                C = beta*C + alpha*A*transpose(B);
            } else if (transA==Trans && transB==Trans) {
                C = beta*C + alpha*transpose(A)*transpose(B);
            } else if (transA==ConjTrans && transB==Trans) {
                C = beta*C + alpha*conjTrans(A)*transpose(B);
            } else if (transA==NoTrans && transB==ConjTrans) {
                C = beta*C + alpha*A*conjTrans(B);
            } else if (transA==Trans && transB==ConjTrans) {
                C = beta*C + alpha*transpose(A)*conjTrans(B);
            } else if (transA==ConjTrans && transB==ConjTrans) {
                C = beta*C + alpha*conjTrans(A)*conjTrans(B);
            }
#       else
            blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}


} // extern "C"
