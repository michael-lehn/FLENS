#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(strsv)(const char      *UPLO,
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

        cblas_strsv(CBLAS_ORDER::CblasColMajor,
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
                BLAS(xerbla)("STRSV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        STrMatrixConstView  A(SFullConstView(*N, *N, _A, *LDA), upLo, diag);
        SDenseVectorView    x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

        // if you only want to test FLENS-BLAS just call
        blas::sv(trans, A, x);
#   endif
}

void
BLAS(dtrsv)(const char      *UPLO,
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

        cblas_dtrsv(CBLAS_ORDER::CblasColMajor,
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
                BLAS(xerbla)("DTRSV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        DTrMatrixConstView  A(DFullConstView(*N, *N, _A, *LDA), upLo, diag);
        DDenseVectorView    x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

        // if you only want to test FLENS-BLAS just call
        blas::sv(trans, A, x);
#   endif
}

void
BLAS(ctrsv)(const char      *UPLO,
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

        cblas_ctrsv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
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
                BLAS(xerbla)("CTRSV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        CTrMatrixConstView  A(CFullConstView(*N, *N, _A, *LDA), upLo, diag);
        CDenseVectorView    x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

        // if you only want to test FLENS-BLAS just call
        blas::sv(trans, A, x);
#   endif
}

void
BLAS(ztrsv)(const char      *UPLO,
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

        cblas_ztrsv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
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
                BLAS(xerbla)("ZTRSV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(_UPLO);
        Transpose    trans = convertTo<Transpose>(_TRANS);
        Diag         diag  = Diag(_DIAG);

        ZTrMatrixConstView  A(ZFullConstView(*N, *N, _A, *LDA), upLo, diag);
        ZDenseVectorView    x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

        // if you only want to test FLENS-BLAS just call
        blas::sv(trans, A, x);

#   endif
}


} // extern "C"
