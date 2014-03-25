#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const float     *_A,
            const INTEGER   *LDA,
            const float     *BETA,
            float           *_C,
            const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS
    
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        
        cblas_ssyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K, *ALPHA,
                    _A, *LDA, *BETA, _C, *LDC);
    
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
            } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
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
    
        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        SGeMatrixConstView  A = SFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _A,
                                               *LDA);

        SSyMatrixView       C(SFullView(*N, *N, _C, *LDC), upLo);

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
            const double    *_A,
            const INTEGER   *LDA,
            const double    *BETA,
            double          *_C,
            const INTEGER   *LDC)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        
        cblas_dsyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K, *ALPHA,
                    _A, *LDA, *BETA, _C, *LDC);
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
            } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
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
    
        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        DGeMatrixConstView  A = DFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _A,
                                               *LDA);

        DSyMatrixView       C(DFullView(*N, *N, _C, *LDC), upLo);

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
            const cfloat    *_A,
            const INTEGER   *LDA,
            const cfloat    *BETA,
            cfloat          *_C,
            const INTEGER   *LDC)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        
        cblas_csyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K,
                    reinterpret_cast<const float *>(ALPHA),
                    reinterpret_cast<const float *>(_A), *LDA,
                    reinterpret_cast<const float *>(BETA),
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
            } else if (_TRANS!='N' && _TRANS!='T') {
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

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        CGeMatrixConstView  A = CFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _A,
                                               *LDA);

        CSyMatrixView       C(CFullView(*N, *N, _C, *LDC), upLo);

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
            const cdouble   *_A,
            const INTEGER   *LDA,
            const cdouble   *BETA,
            cdouble         *_C,
            const INTEGER   *LDC)
{
        
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        
        cblas_zsyrk(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    *N, *K,
                    reinterpret_cast<const double *>(ALPHA),
                    reinterpret_cast<const double *>(_A), *LDA,
                    reinterpret_cast<const double *>(BETA),
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
            } else if (_TRANS!='N' && _TRANS!='T') {
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

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);

        const bool noTrans = (trans==NoTrans || trans==Conj);

        ZGeMatrixConstView  A = ZFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _A,
                                               *LDA);

        ZSyMatrixView       C(ZFullView(*N, *N, _C, *LDC), upLo);

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
