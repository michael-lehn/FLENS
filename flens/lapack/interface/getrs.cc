#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GETRS         sgetrs_
#   define GETRS_REF     sgetrs
#   define GETRS_NAME    "SGETRS"
#elif DOUBLE
#   define GETRS         dgetrs_
#   define GETRS_REF     dgetrs
#   define GETRS_NAME    "DGETRS"
#elif COMPLEX_SINGLE
#   define GETRS         cgetrs_
#   define GETRS_REF     cgetrs
#   define GETRS_NAME    "CGETRS"
#elif COMPLEX_DOUBLE
#   define GETRS         zgetrs_
#   define GETRS_REF     zgetrs
#   define GETRS_NAME    "ZGETRS"
#endif

using namespace flens;

extern "C" {

void
GETRS_REF(char *TRANS, INT *N, INT *NRHS,
          FLOAT *A, INT *LDA,
          INT *IPIV,
          FLOAT *B, INT *LDB,
          INT *INFO);


void
getrsErrorCheck(char *TRANS, INT *N, INT *NRHS,
                FLOAT *A, INT *LDA,
                INT *IPIV,
                FLOAT *B, INT *LDB,
                INT *INFO)
{
    if ((*TRANS!='N') && (*TRANS!='n') && (*TRANS!='T') && (*TRANS!='t')
     && (*TRANS!='C') && (*TRANS!='c') && (*TRANS!='R') && (*TRANS!='r'))
    {
        *INFO = -1;
        return;
    }
    if (*N<0) {
        *INFO = -2;
        return;
    }
    if (*NRHS<0) {
        *INFO = -3;
        return;
    }
    if (*LDA<std::max(INT(1), *N)) {
        *INFO = -5;
        return;
    }
    if (*LDB<std::max(INT(1), *N)) {
        *INFO = -8;
        return;
    }
}


void
GETRS(char *TRANS, INT *N, INT *NRHS,
      FLOAT *A, INT *LDA,
      INT *IPIV,
      FLOAT *B, INT *LDB,
      INT *INFO)
{
    getrsErrorCheck(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GETRS_NAME, INFO, strlen(GETRS_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifndef USE_GETRS_REF
    Transpose trans = NoTrans;

    if ((*TRANS=='N') || (*TRANS=='n')) {
        trans = NoTrans;
    }
    if ((*TRANS=='T') || (*TRANS=='t')) {
        trans = Trans;
    }
    if ((*TRANS=='C') || (*TRANS=='c')) {
        trans = ConjTrans;
    }
    if ((*TRANS=='R') || (*TRANS=='r')) {
        trans = Conj;
    }

    GeMatrix<FSV>        _A = FSV(*N, *N, A, *LDA);
    DenseVector<IAV>     piv = IAV(*N, IPIV, INT(1));
    GeMatrix<FSV>        _B = FSV(*N, *NRHS, B, *LDB);

    lapack::trs(trans, _A, piv, _B);
#   else
    GETRS_REF(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO);
#   endif
}

} // extern "C"

