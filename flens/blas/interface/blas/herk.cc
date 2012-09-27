#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(cherk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const cfloat    *_A,
            const INTEGER   *LDA,
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
    } else if (_TRANS!='N' && _TRANS!='C') {
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
        BLAS(xerbla)("CHERK ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    CGeMatrixConstView  A = CFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    CHeMatrixView       C(CFullView(*N, *N, _C, *LDC), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto beta   = *BETA;


    if (beta==float(1)) {
        if (trans==NoTrans) {
            C += alpha*A*conjTrans(A);
        } else if (trans==ConjTrans) {
            C += alpha*conjTrans(A)*A;
        }
    } else if (beta==float(0)) {
        if (trans==NoTrans) {
            C = alpha*A*conjTrans(A);
        } else if (trans==ConjTrans) {
            C = alpha*conjTrans(A)*A;
        }
    } else {
        if (trans==NoTrans) {
            C = beta*C + alpha*A*conjTrans(A);
        } else if (trans==ConjTrans) {
            C = beta*C + alpha*conjTrans(A)*A;
        }
    }
#   else
    blas::rk(trans, *ALPHA, A, *BETA, C);
#   endif
}

void
BLAS(zherk)(const char      *UPLO,
            const char      *TRANS,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const cdouble   *_A,
            const INTEGER   *LDA,
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
    } else if (_TRANS!='N' && _TRANS!='C') {
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
        BLAS(xerbla)("ZHERK ", &info);
        return;
    }

    StorageUpLo  upLo  = StorageUpLo(_UPLO);
    Transpose    trans = convertTo<Transpose>(_TRANS);

    const bool noTrans = (trans==NoTrans || trans==Conj);

    ZGeMatrixConstView  A = ZFullConstView(noTrans ? *N : *K,
                                           noTrans ? *K : *N,
                                           _A,
                                           *LDA);

    ZHeMatrixView       C(ZFullView(*N, *N, _C, *LDC), upLo);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha  = *ALPHA;
    const auto beta   = *BETA;

    if (beta==double(1)) {
        if (trans==NoTrans) {
            C += alpha*A*conjTrans(A);
        } else if (trans==ConjTrans) {
            C += alpha*conjTrans(A)*A;
        }
    } else if (beta==double(0)) {
        if (trans==NoTrans) {
            C = alpha*A*conjTrans(A);
        } else if (trans==ConjTrans) {
            C = alpha*conjTrans(A)*A;
        }
    } else {
        if (trans==NoTrans) {
            C = beta*C + alpha*A*conjTrans(A);
        } else if (trans==ConjTrans) {
            C = beta*C + alpha*conjTrans(A)*A;
        }
    }
#   else
    blas::rk(trans, *ALPHA, A, *BETA, C);
#   endif
}

} // extern "C"
