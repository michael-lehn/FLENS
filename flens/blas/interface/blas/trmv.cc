#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(strmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const float     *_A,
            const INTEGER   *LDA,
            float           *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_strmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    _A, *LDA, X, *INCX);
    
#   else
    
        using std::abs;
        using std::max;

        char    _UPLO  = toupper(*UPLO);
        char    _TRANS = toupper(*TRANS);
        char    _DIAG  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
                info = 2;
            } else if (_DIAG!='U' && _DIAG!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 6;
            } else if (*INCX==0) {
                info = 8;
            }
            if (info!=0) {
                BLAS(xerbla)("STRMV ", &info);
                return;
            }
#       endif
    
        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        STrMatrixConstView  A(SFullConstView(*N, *N, _A, *LDA), upLo, diag);
        SDenseVectorView    x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            if (trans==NoTrans) {
                x = A*x;
            } else if (trans==Trans || trans==ConjTrans) {
                x = transpose(A)*x;
            }
#       else
            blas::mv(trans, A, x);
#       endif
#   endif
}

void
BLAS(dtrmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const double    *_A,
            const INTEGER   *LDA,
            double          *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS
    
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_dtrmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag), *N,
                    _A, *LDA, X, *INCX);
    
#   else
    
        using std::abs;
        using std::max;

        char    _UPLO  = toupper(*UPLO);
        char    _TRANS = toupper(*TRANS);
        char    _DIAG  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
                info = 2;
            } else if (_DIAG!='U' && _DIAG!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 6;
            } else if (*INCX==0) {
                info = 8;
            }
            if (info!=0) {
                BLAS(xerbla)("DTRMV ", &info);
                return;
            }
#       endif
    
        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        DTrMatrixConstView  A(DFullConstView(*N, *N, _A, *LDA), upLo, diag);
        DDenseVectorView    x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            if (trans==NoTrans) {
                x = A*x;
            } else if (trans==Trans || trans==ConjTrans) {
                x = transpose(A)*x;
            }
#       else
            blas::mv(trans, A, x);
#       endif
#   endif
}

void
BLAS(ctrmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const cfloat    *_A,
            const INTEGER   *LDA,
            cfloat          *X,
            const INTEGER   *INCX)
{
    
#   ifdef TEST_DIRECT_CBLAS
    
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_ctrmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag), *N,
                    reinterpret_cast<const float *>(_A), *LDA,
                    reinterpret_cast<float *>(X), *INCX);
    
#   else
    
    
        using std::abs;
        using std::max;

        char    _UPLO  = toupper(*UPLO);
        char    _TRANS = toupper(*TRANS);
        char    _DIAG  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
                info = 2;
            } else if (_DIAG!='U' && _DIAG!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 6;
            } else if (*INCX==0) {
                info = 8;
            }
            if (info!=0) {
                BLAS(xerbla)("CTRMV ", &info);
                return;
            }
#       endif
    
        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        CTrMatrixConstView  A(CFullConstView(*N, *N, _A, *LDA), upLo, diag);
        CDenseVectorView    x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            if (trans==NoTrans) {
                x = A*x;
            } else if (trans==Trans) {
                x = transpose(A)*x;
            } else if (trans==ConjTrans) {
                x = conjTrans(A)*x;
            }
#       else
            blas::mv(trans, A, x);
#       endif
#   endif
}

void
BLAS(ztrmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const cdouble   *_A,
            const INTEGER   *LDA,
            cdouble         *X,
            const INTEGER   *INCX)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans  = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_ztrmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag), *N,
                    reinterpret_cast<const double *>(_A), *LDA,
                    reinterpret_cast<double *>(X), *INCX);
    
#   else
    
        using std::abs;
        using std::max;

        char    _UPLO  = toupper(*UPLO);
        char    _TRANS = toupper(*TRANS);
        char    _DIAG  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (_UPLO!='U' && _UPLO!='L') {
                info = 1;
            } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
                info = 2;
            } else if (_DIAG!='U' && _DIAG!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*LDA<max(INTEGER(1),*N)) {
                info = 6;
            } else if (*INCX==0) {
                info = 8;
            }
            if (info!=0) {
                BLAS(xerbla)("ZTRMV ", &info);
                return;
            }
#       endif
    
        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        ZTrMatrixConstView  A(ZFullConstView(*N, *N, _A, *LDA), upLo, diag);
        ZDenseVectorView    x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

#       ifdef TEST_OVERLOADED_OPERATORS
            if (trans==NoTrans) {
                x = A*x;
            } else if (trans==Trans) {
                x = transpose(A)*x;
            } else if (trans==ConjTrans) {
                x = conjTrans(A)*x;
            }
#       else
            blas::mv(trans, A, x);
#       endif
#   endif
}


} // extern "C"
