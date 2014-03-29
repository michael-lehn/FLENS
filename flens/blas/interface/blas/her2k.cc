#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(cher2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const cfloat    *ALPHA,
             const cfloat    *_A,
             const INTEGER   *LDA,
             const cfloat    *_B,
             const INTEGER   *LDB,
             const float     *BETA,
             cfloat          *_C,
             const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS

        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);

        cblas_cher2k(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(_A), *LDA,
                    reinterpret_cast<const float *>(_B), *LDB,
                    *BETA,
                    reinterpret_cast<float *>(_C), *LDC);

#   else

        using std::abs;
        using std::max;


        char    _UPLO  = toupper(*UPLO);
        char    _TRANS = toupper(*TRANS);
        INTEGER nRowA  = (_TRANS=='N') ? *N : *K;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (_TRANS!='N' && _TRANS!='C') {
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

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        CGeMatrixConstView  A = CFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _A,
                                               *LDA);

        CGeMatrixConstView  B = CFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _B,
                                               *LDB);

        CHeMatrixView       C(CFullView(*N, *N, _C, *LDC), upLo);

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
             const cdouble   *_A,
             const INTEGER   *LDA,
             const cdouble   *_B,
             const INTEGER   *LDB,
             const double    *BETA,
             cdouble         *_C,
             const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS

        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);

        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);

        cblas_zher2k(CBLAS_ORDER::CblasColMajor,
                     cxxblas::CBLAS::getCblasType(upLo),
                     cxxblas::CBLAS::getCblasType(trans),
                     *N, *K,
                     reinterpret_cast<const double *>(ALPHA),
                     reinterpret_cast<const double *>(_A), *LDA,
                     reinterpret_cast<const double *>(_B), *LDB,
                     *BETA,
                     reinterpret_cast<double *>(_C), *LDC);

#   else

        using std::abs;
        using std::max;

        char    _UPLO  = toupper(*UPLO);
        char    _TRANS = toupper(*TRANS);
        INTEGER nRowA  = (_TRANS=='N') ? *N : *K;

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (_TRANS!='N' && _TRANS!='C') {
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

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        ZGeMatrixConstView  A = ZFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _A,
                                               *LDA);

        ZGeMatrixConstView  B = ZFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _B,
                                               *LDB);

        ZHeMatrixView       C(ZFullView(*N, *N, _C, *LDC), upLo);

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
