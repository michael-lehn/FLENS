#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(stpsv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const float     *AP,
            float           *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_stpsv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    AP,
                    X, *INCX);

#   else
        using std::abs;
        using std::max;

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        char    DIAG_  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
                info = 2;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("STPSV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        STpMatrixConstView    A(SPackedConstView(*N, AP), upLo, diag);
        SDenseVectorView      x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::sv(trans, A, x);
#   endif
}

void
BLAS(dtpsv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const double    *AP,
            double          *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

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

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        char    DIAG_  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
                info = 2;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("DTPSV ", &info);
                return;
            }
#       endif
        StorageUpLo  upLo = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        DTpMatrixConstView    A(DPackedConstView(*N, AP), upLo, diag);
        DDenseVectorView      x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::sv(trans, A, x);
#   endif
}

void
BLAS(ctpsv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const cfloat    *AP,
            cfloat          *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

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

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        char    DIAG_  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
                info = 2;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("CTPSV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        CTpMatrixConstView    A(CPackedConstView(*N, AP), upLo, diag);
        CDenseVectorView      x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::sv(trans, A, x);
#   endif
}

void
BLAS(ztpsv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const cdouble   *AP,
            cdouble         *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

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

        char    UPLO_  = toupper(*UPLO);
        char    TRANS_ = toupper(*TRANS);
        char    DIAG_  = toupper(*DIAG);

#       ifndef NO_INPUT_CHECK
            INTEGER info  = 0;
            if (UPLO_!='U' && UPLO_!='L') {
                info = 1;
            } else if (TRANS_!='N' && TRANS_!='T' && TRANS_!='C') {
                info = 2;
            } else if (DIAG_!='U' && DIAG_!='N') {
                info = 3;
            } else if (*N<0) {
                info = 4;
            } else if (*INCX==0) {
                info = 7;
            }
            if (info!=0) {
                BLAS(xerbla)("ZTPSV ", &info);
                return;
            }
#       endif
        StorageUpLo  upLo = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        ZTpMatrixConstView    A(ZPackedConstView(*N, AP), upLo, diag);
        ZDenseVectorView      x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

        blas::sv(trans, A, x);
#   endif
}

} // extern "C"
