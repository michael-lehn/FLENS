#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyr2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const float     *ALPHA,
             const float     *_A,
             const INTEGER   *LDA,
             const float     *_B,
             const INTEGER   *LDB,
             const float     *BETA,
             float           *_C,
             const INTEGER   *LDC)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
    
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
    
        cblas_ssyr2k(CBLAS_ORDER::CblasColMajor,
                     cxxblas::CBLAS::getCblasType(upLo),
                     cxxblas::CBLAS::getCblasType(trans),
                     *N, *K, *ALPHA,
                     _A, *LDA, _B, *LDB, *BETA, _C, *LDC);
        
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
            } else if (*LDB<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("SSYR2K", &info);
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

        SGeMatrixConstView  B = SFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _B,
                                               *LDB);

        SSyMatrixView       C(SFullView(*N, *N, _C, *LDC), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==float(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==ConjTrans || trans==Trans) {
                    C += alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else if (beta==float(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==ConjTrans || trans==Trans) {
                    C = alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==ConjTrans || trans==Trans) {
                    C = beta*C + alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            }
#       else
            blas::r2k(trans, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}


void
BLAS(dsyr2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const double    *ALPHA,
             const double    *_A,
             const INTEGER   *LDA,
             const double    *_B,
             const INTEGER   *LDB,
             const double    *BETA,
             double          *_C,
             const INTEGER   *LDC)
{
        
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        
        cblas_dsyr2k(CBLAS_ORDER::CblasColMajor,
                     cxxblas::CBLAS::getCblasType(upLo),
                     cxxblas::CBLAS::getCblasType(trans),
                     *N, *K, *ALPHA,
                     _A, *LDA, _B, *LDB, *BETA, _C, *LDC);
        
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
            } else if (*LDB<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("DSYR2K", &info);
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

        DGeMatrixConstView  B = DFullConstView(noTrans ? *N : *K,
                                               noTrans ? *K : *N,
                                               _B,
                                               *LDB);

        DSyMatrixView       C(DFullView(*N, *N, _C, *LDC), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==float(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==ConjTrans || trans==Trans) {
                    C += alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else if (beta==float(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==ConjTrans || trans==Trans) {
                    C = alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==ConjTrans || trans==Trans) {
                    C = beta*C + alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            }
#       else
            blas::r2k(trans, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}

void
BLAS(csyr2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const cfloat    *ALPHA,
             const cfloat    *_A,
             const INTEGER   *LDA,
             const cfloat    *_B,
             const INTEGER   *LDB,
             const cfloat    *BETA,
             cfloat          *_C,
             const INTEGER   *LDC)
{
        
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        
        cblas_csyr2k(CBLAS_ORDER::CblasColMajor,
                     cxxblas::CBLAS::getCblasType(upLo),
                     cxxblas::CBLAS::getCblasType(trans),
                     *N, *K,
                     reinterpret_cast<const float *>(ALPHA),
                     reinterpret_cast<const float *>(_A), *LDA,
                     reinterpret_cast<const float *>(_B), *LDB,
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
            } else if (*LDB<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("CSYR2K", &info);
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

        CSyMatrixView       C(CFullView(*N, *N, _C, *LDC), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==cfloat(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==Trans) {
                    C += alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else if (beta==cfloat(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==Trans) {
                    C = alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==Trans) {
                    C = beta*C + alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            }
#       else
            blas::r2k(trans, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}

void
BLAS(zsyr2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const cdouble   *ALPHA,
             const cdouble   *_A,
             const INTEGER   *LDA,
             const cdouble   *_B,
             const INTEGER   *LDB,
             const cdouble   *BETA,
             cdouble         *_C,
             const INTEGER   *LDC)
{
#   ifdef TEST_DIRECT_CBLAS
    
    char    _UPLO   = toupper(*UPLO);
    char    _TRANS  = toupper(*TRANS);
    
    StorageUpLo    upLo   = StorageUpLo(_UPLO);
    Transpose      trans  = convertTo<Transpose>(_TRANS);
    
    cblas_zsyr2k(CBLAS_ORDER::CblasColMajor,
                 cxxblas::CBLAS::getCblasType(upLo),
                 cxxblas::CBLAS::getCblasType(trans),
                 *N, *K,
                 reinterpret_cast<const double *>(ALPHA),
                 reinterpret_cast<const double *>(_A), *LDA,
                 reinterpret_cast<const double *>(_B), *LDB,
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
            } else if (*LDB<max(INTEGER(1),nRowA)) {
                info = 9;
            } else if (*LDC<max(INTEGER(1),*N)) {
                info = 12;
            }
            if (info!=0) {
                BLAS(xerbla)("ZSYR2K", &info);
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

        ZSyMatrixView       C(ZFullView(*N, *N, _C, *LDC), upLo);

#       ifdef TEST_OVERLOADED_OPERATORS
            const auto alpha  = *ALPHA;
            const auto beta   = *BETA;

            if (beta==cdouble(1)) {
                if (trans==NoTrans) {
                    C += alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==Trans) {
                    C += alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else if (beta==cdouble(0)) {
                if (trans==NoTrans) {
                    C = alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==Trans) {
                    C = alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            } else {
                if (trans==NoTrans) {
                    C = beta*C + alpha*A*transpose(B) + alpha*B*transpose(A);
                } else if (trans==Trans) {
                    C = beta*C + alpha*transpose(A)*B + alpha*transpose(B)*A;
                }
            }
#       else
            blas::r2k(trans, *ALPHA, A, B, *BETA, C);
#       endif
#   endif
}

} // extern "C"
