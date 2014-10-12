#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(strsv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const float     *A_,
            const INTEGER   *LDA,
            float           *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_strsv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    A_, *LDA, X, *INCX);

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

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        STrMatrixConstView  A(SFullConstView(*N, *N, A_, *LDA), upLo, diag);
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
            const double    *A_,
            const INTEGER   *LDA,
            double          *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_dtrsv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    A_, *LDA, X, *INCX);

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

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        DTrMatrixConstView  A(DFullConstView(*N, *N, A_, *LDA), upLo, diag);
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
            const cfloat    *A_,
            const INTEGER   *LDA,
            cfloat          *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_ctrsv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    reinterpret_cast<const float *>(A_), *LDA,
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

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        CTrMatrixConstView  A(CFullConstView(*N, *N, A_, *LDA), upLo, diag);
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
            const cdouble   *A_,
            const INTEGER   *LDA,
            cdouble         *X,
            const INTEGER   *INCX)
{
#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_ztrsv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N,
                    reinterpret_cast<const double *>(A_), *LDA,
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

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        ZTrMatrixConstView  A(ZFullConstView(*N, *N, A_, *LDA), upLo, diag);
        ZDenseVectorView    x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

        // if you only want to test FLENS-BLAS just call
        blas::sv(trans, A, x);

#   endif
}


} // extern "C"
