#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(sgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const float     *ALPHA,
            const float     *_A,
            const INTEGER   *LDA,
            const float     *_B,
            const INTEGER   *LDB,
            const float     *BETA,
            float           *_C,
            const INTEGER   *LDC)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _TRANSA = toupper(*TRANSA);
    char    _TRANSB = toupper(*TRANSB);
    INTEGER nRowA   = (_TRANSA=='N') ? *M : *K;
    INTEGER nRowB   = (_TRANSB=='N') ? *K : *N;

    if (_TRANSA!='N' && _TRANSA!='C' && _TRANSA!='T') {
        info = 1;
    } else if (_TRANSB!='N' && _TRANSB!='C' && _TRANSB!='T') {
        info = 2;
    } else if (*M<0) {
        info = 3;
    } else if (*N<0) {
        info = 4;
    } else if (*K<0) {
        info = 5;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 8;
    } else if (*LDB<max(INTEGER(1),nRowB)) {
        info = 10;
    } else if (*LDC<max(INTEGER(1),*M)) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("SGEMM ", &info);
        return;
    }

    Transpose  transA = convertTo<Transpose>(_TRANSA);
    Transpose  transB = convertTo<Transpose>(_TRANSB);

    const bool noTransA = (transA==NoTrans || transA==Conj);
    const bool noTransB = (transB==NoTrans || transB==Conj);

    const SGeMatrixConstView  A = SFullConstView(noTransA ? *M : *K,
                                                 noTransA ? *K : *M,
                                                 _A,
                                                 *LDA);

    const SGeMatrixConstView  B = SFullConstView(noTransB ? *K : *N,
                                                 noTransB ? *N : *K,
                                                 _B,
                                                 *LDB);

    SGeMatrixView       C = SFullView(*M, *N, _C, *LDC);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;
    const auto beta  = *BETA;

    if (transA==NoTrans && transB==NoTrans) {
        C = beta*C + alpha*A*B;
    } else if (transA==Trans && transB==NoTrans) {
        C = beta*C + alpha*transpose(A)*B;
    } else if (transA==ConjTrans && transB==NoTrans) {
        C = beta*C + alpha*conjTrans(A)*B;
    } else if (transA==NoTrans && transB==Trans) {
        C = beta*C + alpha*A*transpose(B);
    } else if (transA==Trans && transB==Trans) {
        C = beta*C + alpha*transpose(A)*transpose(B);
    } else if (transA==ConjTrans && transB==Trans) {
        C = beta*C + alpha*conjTrans(A)*transpose(B);
    } else if (transA==NoTrans && transB==ConjTrans) {
        C = beta*C + alpha*A*conjTrans(B);
    } else if (transA==Trans && transB==ConjTrans) {
        C = beta*C + alpha*transpose(A)*conjTrans(B);
    } else if (transA==ConjTrans && transB==ConjTrans) {
        C = beta*C + alpha*conjTrans(A)*conjTrans(B);
    }
#   else
    blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#   endif
}


void
BLAS(dgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const double    *ALPHA,
            const double    *_A,
            const INTEGER   *LDA,
            const double    *_B,
            const INTEGER   *LDB,
            const double    *BETA,
            double          *_C,
            const INTEGER   *LDC)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _TRANSA = toupper(*TRANSA);
    char    _TRANSB = toupper(*TRANSB);
    INTEGER nRowA   = (_TRANSA=='N') ? *M : *K;
    INTEGER nRowB   = (_TRANSB=='N') ? *K : *N;

    if (_TRANSA!='N' && _TRANSA!='C' && _TRANSA!='T') {
        info = 1;
    } else if (_TRANSB!='N' && _TRANSB!='C' && _TRANSB!='T') {
        info = 2;
    } else if (*M<0) {
        info = 3;
    } else if (*N<0) {
        info = 4;
    } else if (*K<0) {
        info = 5;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 8;
    } else if (*LDB<max(INTEGER(1),nRowB)) {
        info = 10;
    } else if (*LDC<max(INTEGER(1),*M)) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("DGEMM ", &info);
        return;
    }



    Transpose  transA = convertTo<Transpose>(_TRANSA);
    Transpose  transB = convertTo<Transpose>(_TRANSB);

    const bool noTransA = (transA==NoTrans || transA==Conj);
    const bool noTransB = (transB==NoTrans || transB==Conj);

    const DGeMatrixConstView  A = DFullConstView(noTransA ? *M : *K,
                                                 noTransA ? *K : *M,
                                                 _A,
                                                 *LDA);

    const DGeMatrixConstView  B = DFullConstView(noTransB ? *K : *N,
                                                 noTransB ? *N : *K,
                                                 _B,
                                                 *LDB);

    DGeMatrixView  C = DFullView(*M, *N, _C, *LDC);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;
    const auto beta  = *BETA;

    if (transA==NoTrans && transB==NoTrans) {
        C = beta*C + alpha*A*B;
    } else if (transA==Trans && transB==NoTrans) {
        C = beta*C + alpha*transpose(A)*B;
    } else if (transA==ConjTrans && transB==NoTrans) {
        C = beta*C + alpha*conjTrans(A)*B;
    } else if (transA==NoTrans && transB==Trans) {
        C = beta*C + alpha*A*transpose(B);
    } else if (transA==Trans && transB==Trans) {
        C = beta*C + alpha*transpose(A)*transpose(B);
    } else if (transA==ConjTrans && transB==Trans) {
        C = beta*C + alpha*conjTrans(A)*transpose(B);
    } else if (transA==NoTrans && transB==ConjTrans) {
        C = beta*C + alpha*A*conjTrans(B);
    } else if (transA==Trans && transB==ConjTrans) {
        C = beta*C + alpha*transpose(A)*conjTrans(B);
    } else if (transA==ConjTrans && transB==ConjTrans) {
        C = beta*C + alpha*conjTrans(A)*conjTrans(B);
    }
#   else
    blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#   endif
}

void
BLAS(cgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const cfloat    *ALPHA,
            const cfloat    *_A,
            const INTEGER   *LDA,
            const cfloat    *_B,
            const INTEGER   *LDB,
            const cfloat    *BETA,
            cfloat          *_C,
            const INTEGER   *LDC)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _TRANSA = toupper(*TRANSA);
    char    _TRANSB = toupper(*TRANSB);
    INTEGER nRowA   = (_TRANSA=='N') ? *M : *K;
    INTEGER nRowB   = (_TRANSB=='N') ? *K : *N;

    if (_TRANSA!='N' && _TRANSA!='C' && _TRANSA!='T') {
        info = 1;
    } else if (_TRANSB!='N' && _TRANSB!='C' && _TRANSB!='T') {
        info = 2;
    } else if (*M<0) {
        info = 3;
    } else if (*N<0) {
        info = 4;
    } else if (*K<0) {
        info = 5;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 8;
    } else if (*LDB<max(INTEGER(1),nRowB)) {
        info = 10;
    } else if (*LDC<max(INTEGER(1),*M)) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("CGEMM ", &info);
        return;
    }

    Transpose  transA = convertTo<Transpose>(_TRANSA);
    Transpose  transB = convertTo<Transpose>(_TRANSB);

    const bool noTransA = (transA==NoTrans || transA==Conj);
    const bool noTransB = (transB==NoTrans || transB==Conj);

    const CGeMatrixConstView  A = CFullConstView(noTransA ? *M : *K,
                                                 noTransA ? *K : *M,
                                                 _A,
                                                 *LDA);

    const CGeMatrixConstView  B = CFullConstView(noTransB ? *K : *N,
                                                 noTransB ? *N : *K,
                                                 _B,
                                                 *LDB);

    CGeMatrixView       C = CFullView(*M, *N, _C, *LDC);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;
    const auto beta  = *BETA;

    if (transA==NoTrans && transB==NoTrans) {
        C = beta*C + alpha*A*B;
    } else if (transA==Trans && transB==NoTrans) {
        C = beta*C + alpha*transpose(A)*B;
    } else if (transA==ConjTrans && transB==NoTrans) {
        C = beta*C + alpha*conjTrans(A)*B;
    } else if (transA==NoTrans && transB==Trans) {
        C = beta*C + alpha*A*transpose(B);
    } else if (transA==Trans && transB==Trans) {
        C = beta*C + alpha*transpose(A)*transpose(B);
    } else if (transA==ConjTrans && transB==Trans) {
        C = beta*C + alpha*conjTrans(A)*transpose(B);
    } else if (transA==NoTrans && transB==ConjTrans) {
        C = beta*C + alpha*A*conjTrans(B);
    } else if (transA==Trans && transB==ConjTrans) {
        C = beta*C + alpha*transpose(A)*conjTrans(B);
    } else if (transA==ConjTrans && transB==ConjTrans) {
        C = beta*C + alpha*conjTrans(A)*conjTrans(B);
    }
#   else
    blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#   endif
}


void
BLAS(zgemm)(const char      *TRANSA,
            const char      *TRANSB,
            const INTEGER   *M,
            const INTEGER   *N,
            const INTEGER   *K,
            const cdouble   *ALPHA,
            const cdouble   *_A,
            const INTEGER   *LDA,
            const cdouble   *_B,
            const INTEGER   *LDB,
            const cdouble   *BETA,
            cdouble         *_C,
            const INTEGER   *LDC)
{
    using std::abs;
    using std::max;

    INTEGER info    = 0;
    char    _TRANSA = toupper(*TRANSA);
    char    _TRANSB = toupper(*TRANSB);
    INTEGER nRowA   = (_TRANSA=='N') ? *M : *K;
    INTEGER nRowB   = (_TRANSB=='N') ? *K : *N;

    if (_TRANSA!='N' && _TRANSA!='C' && _TRANSA!='T') {
        info = 1;
    } else if (_TRANSB!='N' && _TRANSB!='C' && _TRANSB!='T') {
        info = 2;
    } else if (*M<0) {
        info = 3;
    } else if (*N<0) {
        info = 4;
    } else if (*K<0) {
        info = 5;
    } else if (*LDA<max(INTEGER(1),nRowA)) {
        info = 8;
    } else if (*LDB<max(INTEGER(1),nRowB)) {
        info = 10;
    } else if (*LDC<max(INTEGER(1),*M)) {
        info = 13;
    }
    if (info!=0) {
        BLAS(xerbla)("ZGEMM ", &info);
        return;
    }



    Transpose  transA = convertTo<Transpose>(_TRANSA);
    Transpose  transB = convertTo<Transpose>(_TRANSB);

    const bool noTransA = (transA==NoTrans || transA==Conj);
    const bool noTransB = (transB==NoTrans || transB==Conj);

    const ZGeMatrixConstView  A = ZFullConstView(noTransA ? *M : *K,
                                                 noTransA ? *K : *M,
                                                 _A,
                                                 *LDA);

    const ZGeMatrixConstView  B = ZFullConstView(noTransB ? *K : *N,
                                                 noTransB ? *N : *K,
                                                 _B,
                                                 *LDB);

    ZGeMatrixView  C = ZFullView(*M, *N, _C, *LDC);

#   ifdef TEST_OVERLOADED_OPERATORS
    const auto alpha = *ALPHA;
    const auto beta  = *BETA;

    if (transA==NoTrans && transB==NoTrans) {
        C = beta*C + alpha*A*B;
    } else if (transA==Trans && transB==NoTrans) {
        C = beta*C + alpha*transpose(A)*B;
    } else if (transA==ConjTrans && transB==NoTrans) {
        C = beta*C + alpha*conjTrans(A)*B;
    } else if (transA==NoTrans && transB==Trans) {
        C = beta*C + alpha*A*transpose(B);
    } else if (transA==Trans && transB==Trans) {
        C = beta*C + alpha*transpose(A)*transpose(B);
    } else if (transA==ConjTrans && transB==Trans) {
        C = beta*C + alpha*conjTrans(A)*transpose(B);
    } else if (transA==NoTrans && transB==ConjTrans) {
        C = beta*C + alpha*A*conjTrans(B);
    } else if (transA==Trans && transB==ConjTrans) {
        C = beta*C + alpha*transpose(A)*conjTrans(B);
    } else if (transA==ConjTrans && transB==ConjTrans) {
        C = beta*C + alpha*conjTrans(A)*conjTrans(B);
    }
#   else
    blas::mm(transA, transB, *ALPHA, A, B, *BETA, C);
#   endif
}


} // extern "C"
