#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(stbsv)(const char      *UPLO,
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
        BLAS(xerbla)("STBSV ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tbsv(ColMajor, upLo,
                  trans, diag,
                  *N, *K,
                  _A, *LDA,
                  X, *INCX);
}

void
BLAS(dtbsv)(const char      *UPLO,
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
        BLAS(xerbla)("DTBSV ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tbsv(ColMajor, upLo,
                  trans, diag,
                  *N, *K,
                  _A, *LDA,
                  X, *INCX);
}

void
BLAS(ctbsv)(const char      *UPLO,
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
        BLAS(xerbla)("CTBSV ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tbsv(ColMajor, upLo,
                  trans, diag,
                  *N, *K,
                  _A, *LDA,
                  X, *INCX);
}

void
BLAS(ztbsv)(const char      *UPLO,
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
        BLAS(xerbla)("ZTBSV ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);
    Diag         diag  = Diag(_DIAG);

    // Until we have band matrices back in FLENS we jsut call cxxblas
    // directly.

    cxxblas::tbsv(ColMajor, upLo,
                  trans, diag,
                  *N, *K,
                  _A, *LDA,
                  X, *INCX);
}


} // extern "C"
