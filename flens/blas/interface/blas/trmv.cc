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
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);
    char    _TRANS = toupper(*TRANS);
    char    _DIAG  = toupper(*DIAG);

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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    STrMatrixConstView  A(SFullConstView(*N, *N, _A, *LDA), upLo, diag);
    SDenseVectorView    x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    if (trans==NoTrans) {
        x = A*x;
    } else if (trans==Trans || trans==ConjTrans) {
        x = transpose(A)*x;
    }
#   else
    blas::mv(trans, A, x);
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
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);
    char    _TRANS = toupper(*TRANS);
    char    _DIAG  = toupper(*DIAG);

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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    DTrMatrixConstView  A(DFullConstView(*N, *N, _A, *LDA), upLo, diag);
    DDenseVectorView    x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    if (trans==NoTrans) {
        x = A*x;
    } else if (trans==Trans || trans==ConjTrans) {
        x = transpose(A)*x;
    }
#   else
    blas::mv(trans, A, x);
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
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);
    char    _TRANS = toupper(*TRANS);
    char    _DIAG  = toupper(*DIAG);

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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    CTrMatrixConstView  A(CFullConstView(*N, *N, _A, *LDA), upLo, diag);
    CDenseVectorView    x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    if (trans==NoTrans) {
        x = A*x;
    } else if (trans==Trans) {
        x = transpose(A)*x;
    } else if (trans==ConjTrans) {
        x = conjTrans(A)*x;
    }
#   else
    blas::mv(trans, A, x);
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
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);
    char    _TRANS = toupper(*TRANS);
    char    _DIAG  = toupper(*DIAG);

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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    ZTrMatrixConstView  A(ZFullConstView(*N, *N, _A, *LDA), upLo, diag);
    ZDenseVectorView    x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

#   ifdef TEST_OVERLOADED_OPERATORS
    if (trans==NoTrans) {
        x = A*x;
    } else if (trans==Trans) {
        x = transpose(A)*x;
    } else if (trans==ConjTrans) {
        x = conjTrans(A)*x;
    }
#   else
    blas::mv(trans, A, x);
#   endif
}


} // extern "C"
