#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define ORG2R         sorg2r_
#   define ORG2R_REF     sorg2r
#   define ORG2R_NAME    "SORG2R"
#elif DOUBLE
#   define ORG2R         dorg2r_
#   define ORG2R_REF     dorg2r
#   define ORG2R_NAME    "DORG2R"
#elif COMPLEX_SINGLE
#   define ORG2R         corg2r_
#   define ORG2R_REF     corg2r
#   define ORG2R_NAME    "CORG2R"
#elif COMPLEX_DOUBLE
#   define ORG2R         zorg2r_
#   define ORG2R_REF     zorg2r
#   define ORG2R_NAME    "ZORG2R"
#endif

using namespace flens;

extern "C" {

void
ORG2R_REF(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
          FLOAT *WORK, INT *INFO);

void
org2rErrorCheck(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
                FLOAT *WORK, INT *INFO)
{
    if (*M<0) {
        *INFO = -1;
    } else if ((*N<0) || (*N>*M)) {
        *INFO = -2;
    } else if ((*K<0) || (*K>*N)) {
        *INFO = -3;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -5;
    }
}


void
ORG2R(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
      FLOAT *WORK, INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << ORG2R_NAME << std::endl;
#   endif

    org2rErrorCheck(M, N, K, A, LDA, TAU, WORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(ORG2R_NAME, INFO, strlen(ORG2R_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifdef USE_ORG2R_REF
    // call LAPACK implementation

    ORG2R_REF(M, N, K, A, LDA, TAU, WORK, INFO);

#   elif defined(CMP_W_ORG2R_REF)
    // call LAPACK and FLENS implementation, compare results

    if (*N==0) {
        *INFO = 0;
        return;
    }

    GeMatrix<FSV>        _A = FSV(*M, *N, A, *LDA);
    DenseVector<AV>      _tau = AV(*K, TAU);
    DenseVector<AV>      _work = AV(*LWORK, WORK);


    GeMatrix<FS>         __A = _A;
    DenseVector<AR>      __tau = _tau;
    DenseVector<AR>      __work = _work;


    lapack::org2r(*K, __A, __tau, __work);
    ORG2R_REF(M, N, K, A, LDA, TAU, WORK, INFO);

    if (isDifferent(_A,__A)) {
        std::cerr << "_A = " << _A << std::endl;
        std::cerr << "__A = " << __A << std::endl;
        assert(0);
    }

    if (isDifferent(_tau,__tau)) {
        std::cerr << "_tau = " << _tau << std::endl;
        std::cerr << "__tau = " << __tau << std::endl;
        assert(0);
    }

    if (isDifferent(_work,__work)) {
        std::cerr << "ERROR:" << std::endl;
        for (INT i=_work.firstIndex(); i<=_work.lastIndex(); ++i) {
            std::cerr << "LAPACK: _work(" << i << ") = " << _work(i) << std::endl
                      << "FLENS: __work(" << i << ") = " << __work(i)
                      << std::endl << std::endl;
        }
        assert(0);
    }

#   else
    // call FLENS implementation

    GeMatrix<FSV>        _A = FSV(*M, *N, A, *LDA);
    DenseVector<AV>      _tau = AV(*K, TAU);
    DenseVector<AV>      _work = AV(*LWORK, WORK);

    lapack::org2r(*K, _A, _tau, _work);
#   endif
}

} // extern "C"

