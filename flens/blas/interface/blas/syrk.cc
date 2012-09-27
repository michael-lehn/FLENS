#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(ssyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const float     *_A,
            const INTEGER   *LDA,
            const float     *BETA,
            float           *_C,
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
    } else if (*LDC<max(INTEGER(1),*N)) {
        info = 10;
    }
    if (info!=0) {
        BLAS(xerbla)("SSYRK ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    SGeMatrixConstView  A = SFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    SSyMatrixView       C(SFullView(*N, *N, _C, *LDC), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto beta   = *BETA;

    if (beta==float(1)) {
        if (trans==NoTrans) {
            C += alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C += alpha*transpose(A)*A;
        }
    } else if (beta==float(0)) {
        if (trans==NoTrans) {
            C = alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = alpha*transpose(A)*A;
        }
    } else {
        if (trans==NoTrans) {
            C = beta*C + alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = beta*C + alpha*transpose(A)*A;
        }
    }
#   else
    blas::rk(trans, *ALPHA, A, *BETA, C);
#   endif
}

void
BLAS(dsyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const double    *_A,
            const INTEGER   *LDA,
            const double    *BETA,
            double          *_C,
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
    } else if (*LDC<max(INTEGER(1),*N)) {
        info = 10;
    }
    if (info!=0) {
        BLAS(xerbla)("DSYRK ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    DGeMatrixConstView  A = DFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    DSyMatrixView       C(DFullView(*N, *N, _C, *LDC), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto beta   = *BETA;

    if (beta==double(1)) {
        if (trans==NoTrans) {
            C += alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C += alpha*transpose(A)*A;
        }
    } else if (beta==double(0)) {
        if (trans==NoTrans) {
            C = alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = alpha*transpose(A)*A;
        }
    } else {
        if (trans==NoTrans) {
            C = beta*C + alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = beta*C + alpha*transpose(A)*A;
        }
    }
#   else
    blas::rk(trans, *ALPHA, A, *BETA, C);
#   endif
}

void
BLAS(csyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const cfloat    *ALPHA,
            const cfloat    *_A,
            const INTEGER   *LDA,
            const cfloat    *BETA,
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
    } else if (_TRANS!='N' && _TRANS!='T') {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*K<0) {
        info = 4;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 7;
    } else if (*LDC<max(INTEGER(1),*N)) {
        info = 10;
    }
    if (info!=0) {
        BLAS(xerbla)("CSYRK ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    CGeMatrixConstView  A = CFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    CSyMatrixView       C(CFullView(*N, *N, _C, *LDC), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto beta   = *BETA;

    if (beta==cfloat(1)) {
        if (trans==NoTrans) {
            C += alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C += alpha*transpose(A)*A;
        }
    } else if (beta==cfloat(0)) {
        if (trans==NoTrans) {
            C = alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = alpha*transpose(A)*A;
        }
    } else {
        if (trans==NoTrans) {
            C = beta*C + alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = beta*C + alpha*transpose(A)*A;
        }
    }
#   else
    blas::rk(trans, *ALPHA, A, *BETA, C);
#   endif
}

void
BLAS(zsyrk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const cdouble   *ALPHA,
            const cdouble   *_A,
            const INTEGER   *LDA,
            const cdouble   *BETA,
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
    } else if (_TRANS!='N' && _TRANS!='T') {
        info = 2;
    } else if (*N<0) {
        info = 3;
    } else if (*K<0) {
        info = 4;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 7;
    } else if (*LDC<max(INTEGER(1),*N)) {
        info = 10;
    }
    if (info!=0) {
        BLAS(xerbla)("ZSYRK ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    ZGeMatrixConstView  A = ZFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    ZSyMatrixView       C(ZFullView(*N, *N, _C, *LDC), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto beta   = *BETA;

    if (beta==cdouble(1)) {
        if (trans==NoTrans) {
            C += alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C += alpha*transpose(A)*A;
        }
    } else if (beta==cdouble(0)) {
        if (trans==NoTrans) {
            C = alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = alpha*transpose(A)*A;
        }
    } else {
        if (trans==NoTrans) {
            C = beta*C + alpha*A*transpose(A);
        } else if (trans==Trans || trans==ConjTrans) {
            C = beta*C + alpha*transpose(A)*A;
        }
    }
#   else
    blas::rk(trans, *ALPHA, A, *BETA, C);
#   endif
}


} // extern "C"
