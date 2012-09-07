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
        BLAS(xerbla)("STPSV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tpsv(ColMajor, upLo,
                  trans, diag,
                  *N,
                  AP,
                  X, *INCX);
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
        BLAS(xerbla)("DTPSV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tpsv(ColMajor, upLo,
                  trans, diag,
                  *N,
                  AP,
                  X, *INCX);
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
        BLAS(xerbla)("CTPSV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tpsv(ColMajor, upLo,
                  trans, diag,
                  *N,
                  AP,
                  X, *INCX);
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
        BLAS(xerbla)("ZTPSV ", &info);
        return;
    }

    StorageUpLo  upLo = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tpsv(ColMajor, upLo,
                  trans, diag,
                  *N,
                  AP,
                  X, *INCX);
}

} // extern "C"
