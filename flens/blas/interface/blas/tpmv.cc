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
    } else if (*INCX==0) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("STPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    STpMatrixConstView    A(SPackedConstView(*N, AP), upLo, diag);
    SDenseVectorView      x(SArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
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
    } else if (*INCX==0) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("DTPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    DTpMatrixConstView    A(DPackedConstView(*N, AP), upLo, diag);
    DDenseVectorView      x(DArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
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
    } else if (*INCX==0) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("CTPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    CTpMatrixConstView    A(CPackedConstView(*N, AP), upLo, diag);
    CDenseVectorView      x(CArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
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
    } else if (*INCX==0) {
        info = 7;
    }
    if (info!=0) {
        BLAS(xerbla)("ZTPMV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    ZTpMatrixConstView    A(ZPackedConstView(*N, AP), upLo, diag);
    ZDenseVectorView      x(ZArrayView(*N, X, abs(*INCX)), *INCX<0);

    blas::mv(trans, A, x);
}


} // extern "C"
