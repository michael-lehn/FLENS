#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(cher2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const cfloat    *ALPHA,
             const cfloat    *_A,
             const INTEGER   *LDA,
             const cfloat    *_B,
             const INTEGER   *LDB,
             const float     *BETA,
             cfloat          *_C,
             const INTEGER   *LDC)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);
    char    _TRANS = toupper(*TRANS);
    INTEGER nRowA  = (_TRANS=='N') ? *N : *K;

    if (_UPLO!='U' && _UPLO!='L') {
        info = 1;
    } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*K<0) {
        info = 4;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 7;
    } else if (*LDB<max(INTEGER(1),nRowA)) {
        info = 9;
    } else if (*LDC<max(INTEGER(1),*N)) {
        info = 12;
    }
    if (info!=0) {
        BLAS(xerbla)("CHER2K", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    CGeMatrixConstView  A = CFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    CGeMatrixConstView  B = CFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _B,
                                           *LDB);

    CHeMatrixView       C(CFullView(*N, *N, _C, *LDC), upLo);

    // if you only want to test FLENS-BLAS just call
    blas::r2k(trans, *ALPHA, A, B, *BETA, C);
}


void
BLAS(zher2k)(const char      *UPLO,
             const char      *TRANS,
             const INTEGER   *N,
             const INTEGER   *K,
             const cdouble   *ALPHA,
             const cdouble   *_A,
             const INTEGER   *LDA,
             const cdouble   *_B,
             const INTEGER   *LDB,
             const double    *BETA,
             cdouble         *_C,
             const INTEGER   *LDC)
{
    using std::abs;
    using std::max;

    INTEGER info   = 0;
    char    _UPLO  = toupper(*UPLO);
    char    _TRANS = toupper(*TRANS);
    INTEGER nRowA  = (_TRANS=='N') ? *N : *K;

    if (_UPLO!='U' && _UPLO!='L') {
        info = 1;
    } else if (_TRANS!='N' && _TRANS!='T' && _TRANS!='C') {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*K<0) {
        info = 4;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 7;
    } else if (*LDB<max(INTEGER(1),nRowA)) {
        info = 9;
    } else if (*LDC<max(INTEGER(1),*N)) {
        info = 12;
    }
    if (info!=0) {
        BLAS(xerbla)("ZHER2K", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    ZGeMatrixConstView  A = ZFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    ZGeMatrixConstView  B = ZFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _B,
                                           *LDB);

    ZHeMatrixView       C(ZFullView(*N, *N, _C, *LDC), upLo);

    // if you only want to test FLENS-BLAS just call
    blas::r2k(trans, *ALPHA, A, B, *BETA, C);
}

} // extern "C"
