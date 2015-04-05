#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(strmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const float     *A_,
            const INTEGER   *LDA,
            float           *X,
            const INTEGER   *INCX)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: strmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_strmv(CBLAS_ORDER::CblasColMajor,
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
                BLAS(xerbla)("STRMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        STrMatrixConstView  A(SFullConstView(*N, *N, *LDA, A_), upLo, diag);
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
            const double    *A_,
            const INTEGER   *LDA,
            double          *X,
            const INTEGER   *INCX)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: dtrmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_dtrmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag), *N,
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
                BLAS(xerbla)("DTRMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        DTrMatrixConstView  A(DFullConstView(*N, *N, *LDA, A_), upLo, diag);
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
            const cfloat    *A_,
            const INTEGER   *LDA,
            cfloat          *X,
            const INTEGER   *INCX)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: ctrmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_ctrmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag), *N,
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
                BLAS(xerbla)("CTRMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        CTrMatrixConstView  A(CFullConstView(*N, *N, *LDA, A_), upLo, diag);
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
            const cdouble   *A_,
            const INTEGER   *LDA,
            cdouble         *X,
            const INTEGER   *INCX)
{
    BLAS_DEBUG_OUT("BLAS INTERFACE: ztrmv");

#   ifdef TEST_DIRECT_CBLAS

        char    UPLO_   = toupper(*UPLO);
        char    TRANS_  = toupper(*TRANS);
        char    DIAG_   = toupper(*DIAG);

        StorageUpLo    upLo   = StorageUpLo(UPLO_);
        Transpose      trans  = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_ztrmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag), *N,
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
                BLAS(xerbla)("ZTRMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        ZTrMatrixConstView  A(ZFullConstView(*N, *N, *LDA, A_), upLo, diag);
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
