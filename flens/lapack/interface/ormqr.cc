#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define ORMQR         sormqr_
#   define ORMQR_REF     sormqr
#   define ORMQR_NAME    "SORMQR"
#elif DOUBLE
#   define ORMQR         dormqr_
#   define ORMQR_REF     dormqr
#   define ORMQR_NAME    "DORMQR"
#elif COMPLEX_SINGLE
#   define ORMQR         cormqr_
#   define ORMQR_REF     cormqr
#   define ORMQR_NAME    "CORMQR"
#elif COMPLEX_DOUBLE
#   define ORMQR         zormqr_
#   define ORMQR_REF     zormqr
#   define ORMQR_NAME    "ZORMQR"
#endif

using namespace flens;

extern "C" {

void
ORMQR_REF(char *SIDE, char *TRANS,
          INT *M, INT *N, INT *K,
          FLOAT *A, INT *LDA,
          FLOAT *TAU,
          FLOAT *C, INT *LDC,
          FLOAT *WORK, INT *LWORK,
          INT *INFO);

void
ormqrErrorCheck(char *SIDE, char *TRANS,
                INT *M, INT *N, INT *K,
                FLOAT *A, INT *LDA,
                FLOAT *TAU,
                FLOAT *C, INT *LDC,
                FLOAT *WORK, INT *LWORK,
                INT *INFO)
{
    bool lQuery = (*LWORK==-1);

    bool left = (*SIDE=='l') || (*SIDE=='L');
    bool right = (*SIDE=='r') || (*SIDE=='R');

    INT nq, nw;
    if (left) {
        nq = *M;
        nw = *N;
    } else {
        nq = *N;
        nw = *M;
    }

    bool noTrans = (*TRANS=='N') || (*TRANS=='n')
                || (*TRANS=='R') || (*TRANS=='r');
    bool trans   = (*TRANS=='T') || (*TRANS=='t')
                || (*TRANS=='C') || (*TRANS=='c');

    if ((!left) && (!right)) {
        *INFO = -1;
    } else if ((!noTrans) && (!trans)) {
        *INFO = -2;
    } else if (*M<0) {
        *INFO = -3;
    } else if (*N<0) {
        *INFO = -4;
    } else if ((*K<0) || (*K>nq)) {
        *INFO = -5;
    } else if (*LDA<std::max(INT(1), nq)) {
        *INFO = -7;
    } else if (*LDC<std::max(INT(1), *M)) {
        *INFO = -10;
    } else if ((*LWORK<std::max(INT(1), nw)) && (!lQuery)) {
        *INFO = -12;
    }
}


void
ORMQR(char *SIDE, char *TRANS,
      INT *M, INT *N, INT *K,
      FLOAT *A, INT *LDA,
      FLOAT *TAU,
      FLOAT *C, INT *LDC,
      FLOAT *WORK, INT *LWORK,
      INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << ORMQR_NAME << std::endl;
#   endif

    ormqrErrorCheck(SIDE, TRANS, M, N, K, A, LDA,
                    TAU, C, LDC, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(ORMQR_NAME, INFO, strlen(ORMQR_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_ORMQR_REF
    Side _side = ((*SIDE=='l') || (*SIDE=='L')) ? Left : Right;
    Transpose _trans;
    if ((*TRANS=='N') || (*TRANS=='n')) {
        _trans = NoTrans;
    } else if ((*TRANS=='R') || (*TRANS=='r')) {
        _trans = Conj;
    } else if ((*TRANS=='T') || (*TRANS=='t')) {
        _trans = Trans;
    } else if ((*TRANS=='C') || (*TRANS=='c')) {
        _trans = ConjTrans;
    } else {
        ASSERT(0);
    }

    GeMatrix<FSV>        _A = FSV((_side==Left) ? *M : *N, *K, A, *LDA);
    GeMatrix<FSV>        _C = FSV(*M, *N, C, *LDC);
    DenseVector<AV>      _tau = AV(*K, TAU);

    DenseVector<AV>      _work = AV(*LWORK, WORK);

    // for testing worksize query/resizing pass an empty vector:
    // DenseVector<AR>      __work;

    lapack::ormqr(_side, _trans, _A, _tau, _C, _work);

#   else
    ORMQR_REF(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
#   endif
}

} // extern "C"

