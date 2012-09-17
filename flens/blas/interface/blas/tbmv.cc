#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(stbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    STbMatrixConstView  A(SBandConstView(*N, *N, KL, KU, _A, *LDA), upLo, diag);
    SDenseVectorView    x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
}

void
BLAS(dtbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    DTbMatrixConstView  A(DBandConstView(*N, *N, KL, KU, _A, *LDA), upLo, diag);
    DDenseVectorView    x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
}

void
BLAS(ctbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    CTbMatrixConstView  A(CBandConstView(*N, *N, KL, KU, _A, *LDA), upLo, diag);
    CDenseVectorView    x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
}

void
BLAS(ztbmv)(const char      *UPLO,
            const char      *TRANS,
            const char      *DIAG,
            const INTEGER   *N,
            const INTEGER   *K,
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

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    INTEGER      KL   = (upLo==Lower) ? *K : 0;
    INTEGER      KU   = (upLo==Upper) ? *K : 0;

    ZTbMatrixConstView  A(ZBandConstView(*N, *N, KL, KU, _A, *LDA), upLo, diag);
    ZDenseVectorView    x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
}


} // extern "C"
