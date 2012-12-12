#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(stpmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const float     *AP,
            float           *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_stpmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    AP,
                    X, *INCX);
        
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
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("STPMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        STpMatrixConstView    A(SPackedConstView(*N, AP), upLo, diag);
        SDenseVectorView      x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

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
BLAS(dtpmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const double    *AP,
            double          *X,
            const INTEGER   *INCX)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_dtpmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    AP,
                    X, *INCX);
        
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
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("DTPMV ", &info);
                return;
            }
#       endif
    
        StorageUpLo  upLo = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        DTpMatrixConstView    A(DPackedConstView(*N, AP), upLo, diag);
        DDenseVectorView      x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

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
BLAS(ctpmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const cfloat    *AP,
            cfloat          *X,
            const INTEGER   *INCX)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_ctpmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    reinterpret_cast<const float *>(AP),
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
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("CTPMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        CTpMatrixConstView    A(CPackedConstView(*N, AP), upLo, diag);
        CDenseVectorView      x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

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
BLAS(ztpmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const cdouble   *AP,
            cdouble         *X,
            const INTEGER   *INCX)
{
    
#   ifdef TEST_DIRECT_CBLAS
        
        char    _UPLO   = toupper(*UPLO);
        char    _TRANS  = toupper(*TRANS);
        char    _DIAG   = toupper(*DIAG);
        
        StorageUpLo    upLo   = StorageUpLo(_UPLO);
        Transpose      trans = convertTo<Transpose>(_TRANS);
        Diag           diag   = Diag(_DIAG);
        
        cblas_ztpmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    reinterpret_cast<const double *>(AP),
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
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("ZTPMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        ZTpMatrixConstView    A(ZPackedConstView(*N, AP), upLo, diag);
        ZDenseVectorView      x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

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
