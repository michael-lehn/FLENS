#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(stbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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
        Transpose      trans = convertTo<Transpose>(TRANS_);
        Diag           diag   = Diag(DIAG_);

        cblas_stbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N, *K,
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
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<(*K + 1)) {
                info = 7;
            } else if (*INCX==0) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("STBMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        STbMatrixConstView  A(SBandConstView(*N, *N, KL, KU, A_, *LDA),
                              upLo, diag);
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
BLAS(dtbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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
        Transpose      trans  = convertTo<Transpose>(TRANS_ );
        Diag           diag   = Diag(DIAG_);

        cblas_dtbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N, *K,
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
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<(*K + 1)) {
                info = 7;
            } else if (*INCX==0) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("DTBMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        DTbMatrixConstView  A(DBandConstView(*N, *N, KL, KU, A_, *LDA),
                              upLo, diag);
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
BLAS(ctbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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

        cblas_ctbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N, *K,
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
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<(*K + 1)) {
                info = 7;
            } else if (*INCX==0) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("CTBMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        CTbMatrixConstView  A(CBandConstView(*N, *N, KL, KU, A_, *LDA),
                              upLo, diag);
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
BLAS(ztbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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

        cblas_ztbmv(CBLAS_ORDER::CblasColMajor,
                    cxxblas::CBLAS::getCblasType(upLo),
                    cxxblas::CBLAS::getCblasType(trans),
                    cxxblas::CBLAS::getCblasType(diag),
                    *N, *K,
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
            } else if (*K<0) {
                info = 5;
            } else if (*LDA<(*K + 1)) {
                info = 7;
            } else if (*INCX==0) {
                info = 9;
            }
            if (info!=0) {
                BLAS(xerbla)("ZTBMV ", &info);
                return;
            }
#       endif

        StorageUpLo  upLo  = StorageUpLo(UPLO_);
        Transpose    trans = convertTo<Transpose>(TRANS_);
        Diag         diag  = Diag(DIAG_);

        INTEGER      KL   = (upLo==Lower) ? *K : 0;
        INTEGER      KU   = (upLo==Upper) ? *K : 0;

        ZTbMatrixConstView  A(ZBandConstView(*N, *N, KL, KU, A_, *LDA),
                              upLo, diag);
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
