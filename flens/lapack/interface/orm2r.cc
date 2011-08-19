#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define ORM2R         sorm2r_
#   define ORM2R_REF     sorm2r
#   define ORM2R_NAME    "SORM2R"
#elif DOUBLE
#   define ORM2R         dorm2r_
#   define ORM2R_REF     dorm2r
#   define ORM2R_NAME    "DORM2R"
#elif COMPLEX_SINGLE
#   define ORM2R         corm2r_
#   define ORM2R_REF     corm2r
#   define ORM2R_NAME    "CORM2R"
#elif COMPLEX_DOUBLE
#   define ORM2R         zorm2r_
#   define ORM2R_REF     zorm2r
#   define ORM2R_NAME    "ZORM2R"
#endif

using namespace flens;

extern "C" {

void
ORM2R_REF(char *SIDE, char *TRANS,
          INT *M, INT *N, INT *K,
          FLOAT *A, INT *LDA,
          FLOAT *TAU,
          FLOAT *C, INT *LDC,
          FLOAT *WORK,
          INT *INFO);

void
orm2rErrorCheck(char *SIDE, char *TRANS,
                INT *M, INT *N, INT *K,
                FLOAT *A, INT *LDA,
                FLOAT *TAU,
                FLOAT *C, INT *LDC,
                FLOAT *WORK,
                INT *INFO)
{
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
    }
}


void
ORM2R(char *SIDE, char *TRANS,
      INT *M, INT *N, INT *K,
      FLOAT *A, INT *LDA,
      FLOAT *TAU,
      FLOAT *C, INT *LDC,
      FLOAT *WORK,
      INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << ORM2R_NAME << std::endl;
#   endif

    orm2rErrorCheck(SIDE, TRANS, M, N, K, A, LDA,
                    TAU, C, LDC, WORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(ORM2R_NAME, INFO, strlen(ORM2R_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_ORM2R_REF
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

    DenseVector<AV>      _work = AV((_side==Left) ? *N : *M, WORK);

    // for testing worksize query/resizing pass an empty vector:
    // DenseVector<AR>      __work;

    lapack::orm2r(_side, _trans, _A, _tau, _C, _work);

#   else
    std::cerr << "ORM2R_REF" << std::endl;
    ORM2R_REF(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO);
#   endif
}

} // extern "C"

