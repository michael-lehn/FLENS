#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(strsm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const float     *ALPHA,
            const float     *_A,
            const INTEGER   *LDA,
            float           *_B,
            const INTEGER   *LDB)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _SIDE   = toupper(*SIDE);
    char    _UPLO   = toupper(*UPLO);
    char    _TRANSA = toupper(*TRANSA);
    char    _DIAG   = toupper(*DIAG);
    INTEGER nRowA   = (_SIDE=='L') ? *M : *N;

    if (_SIDE!='L' && _SIDE!='R') {
        info = 1;
    } else if (_UPLO!='U' && _UPLO!='L') {
        info = 2;
    } else if (_TRANSA!='N' && _TRANSA!='T' && _TRANSA!='C') {
        info = 3;
    } else if (_DIAG!='U' && _DIAG!='N') {
        info = 4;
    } else if (*M<0) {
        info = 5;
    } else if (*N<0) {
        info = 6;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 9;
    } else if (*LDB<max(INTEGER(1),*M)) {
        info = 11;
    }
    if (info!=0) {
        BLAS(xerbla)("STRSM ", &info);
        return;
    }

    Side           side   = Side(_SIDE);
    StorageUpLo    upLo   = StorageUpLo(_UPLO);
    Transpose      transA = convertTo<Transpose>(_TRANSA);
    Diag           diag   = Diag(_DIAG);
    const INTEGER  ka     = side==Left ? *M : *N;

    STrMatrixConstView  A(SFullConstView(ka, ka, _A, *LDA), upLo, diag);
    SGeMatrixView       B = SFullView(*M, *N, _B, *LDB);

    blas::sm(side, transA, *ALPHA, A, B);
}

void
BLAS(dtrsm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const double    *ALPHA,
            const double    *_A,
            const INTEGER   *LDA,
            double          *_B,
            const INTEGER   *LDB)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _SIDE   = toupper(*SIDE);
    char    _UPLO   = toupper(*UPLO);
    char    _TRANSA = toupper(*TRANSA);
    char    _DIAG   = toupper(*DIAG);
    INTEGER nRowA   = (_SIDE=='L') ? *M : *N;

    if (_SIDE!='L' && _SIDE!='R') {
        info = 1;
    } else if (_UPLO!='U' && _UPLO!='L') {
        info = 2;
    } else if (_TRANSA!='N' && _TRANSA!='T' && _TRANSA!='C') {
        info = 3;
    } else if (_DIAG!='U' && _DIAG!='N') {
        info = 4;
    } else if (*M<0) {
        info = 5;
    } else if (*N<0) {
        info = 6;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 9;
    } else if (*LDB<max(INTEGER(1),*M)) {
        info = 11;
    }
    if (info!=0) {
        BLAS(xerbla)("DTRSM ", &info);
        return;
    }

    Side           side   = Side(_SIDE);
    StorageUpLo    upLo   = StorageUpLo(_UPLO);
    Transpose      transA = convertTo<Transpose>(_TRANSA);
    Diag           diag   = Diag(_DIAG);
    const INTEGER  ka     = side==Left ? *M : *N;

    DTrMatrixConstView  A(DFullConstView(ka, ka, _A, *LDA), upLo, diag);
    DGeMatrixView       B = DFullView(*M, *N, _B, *LDB);

    blas::sm(side, transA, *ALPHA, A, B);
}

void
BLAS(ctrsm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const cfloat    *ALPHA,
            const cfloat    *_A,
            const INTEGER   *LDA,
            cfloat          *_B,
            const INTEGER   *LDB)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _SIDE   = toupper(*SIDE);
    char    _UPLO   = toupper(*UPLO);
    char    _TRANSA = toupper(*TRANSA);
    char    _DIAG   = toupper(*DIAG);
    INTEGER nRowA   = (_SIDE=='L') ? *M : *N;

    if (_SIDE!='L' && _SIDE!='R') {
        info = 1;
    } else if (_UPLO!='U' && _UPLO!='L') {
        info = 2;
    } else if (_TRANSA!='N' && _TRANSA!='T' && _TRANSA!='C') {
        info = 3;
    } else if (_DIAG!='U' && _DIAG!='N') {
        info = 4;
    } else if (*M<0) {
        info = 5;
    } else if (*N<0) {
        info = 6;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 9;
    } else if (*LDB<max(INTEGER(1),*M)) {
        info = 11;
    }
    if (info!=0) {
        BLAS(xerbla)("CTRSM ", &info);
        return;
    }

    Side           side   = Side(_SIDE);
    StorageUpLo    upLo   = StorageUpLo(_UPLO);
    Transpose      transA = convertTo<Transpose>(_TRANSA);
    Diag           diag   = Diag(_DIAG);
    const INTEGER  ka     = side==Left ? *M : *N;

    CTrMatrixConstView  A(CFullConstView(ka, ka, _A, *LDA), upLo, diag);
    CGeMatrixView       B = CFullView(*M, *N, _B, *LDB);

    blas::sm(side, transA, *ALPHA, A, B);
}

void
BLAS(ztrsm)(const char      *SIDE,
            const char      *UPLO,
            const char      *TRANSA,
            const char      *DIAG,
            const INTEGER   *M,
            const INTEGER   *N,
            const cdouble   *ALPHA,
            const cdouble   *_A,
            const INTEGER   *LDA,
            cdouble         *_B,
            const INTEGER   *LDB)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _SIDE   = toupper(*SIDE);
    char    _UPLO   = toupper(*UPLO);
    char    _TRANSA = toupper(*TRANSA);
    char    _DIAG   = toupper(*DIAG);
    INTEGER nRowA   = (_SIDE=='L') ? *M : *N;

    if (_SIDE!='L' && _SIDE!='R') {
        info = 1;
    } else if (_UPLO!='U' && _UPLO!='L') {
        info = 2;
    } else if (_TRANSA!='N' && _TRANSA!='T' && _TRANSA!='C') {
        info = 3;
    } else if (_DIAG!='U' && _DIAG!='N') {
        info = 4;
    } else if (*M<0) {
        info = 5;
    } else if (*N<0) {
        info = 6;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 9;
    } else if (*LDB<max(INTEGER(1),*M)) {
        info = 11;
    }
    if (info!=0) {
        BLAS(xerbla)("ZTRSM ", &info);
        return;
    }

    Side           side   = Side(_SIDE);
    StorageUpLo    upLo   = StorageUpLo(_UPLO);
    Transpose      transA = convertTo<Transpose>(_TRANSA);
    Diag           diag   = Diag(_DIAG);
    const INTEGER  ka     = side==Left ? *M : *N;

    ZTrMatrixConstView  A(ZFullConstView(ka, ka, _A, *LDA), upLo, diag);
    ZGeMatrixView       B = ZFullView(*M, *N, _B, *LDB);

    blas::sm(side, transA, *ALPHA, A, B);
}


} // extern "C"
