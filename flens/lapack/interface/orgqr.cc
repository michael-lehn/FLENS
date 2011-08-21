#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define ORGQR         sorgqr_
#   define ORGQR_REF     sorgqr
#   define ORGQR_NAME    "SORGQR"
#elif DOUBLE
#   define ORGQR         dorgqr_
#   define ORGQR_REF     dorgqr
#   define ORGQR_NAME    "DORGQR"
#elif COMPLEX_SINGLE
#   define ORGQR         corgqr_
#   define ORGQR_REF     corgqr
#   define ORGQR_NAME    "CORGQR"
#elif COMPLEX_DOUBLE
#   define ORGQR         zorgqr_
#   define ORGQR_REF     zorgqr
#   define ORGQR_NAME    "ZORGQR"
#endif

using namespace flens;

extern "C" {

void
ORGQR_REF(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
          FLOAT *WORK, INT *LWORK, INT *INFO);

void
orgqrErrorCheck(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
                FLOAT *WORK, INT *LWORK, INT *INFO)
{
    bool lQuery = (*LWORK==-1);
    
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*K<0) {
        *INFO = -3;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -5;
    } else if ((*LWORK<std::max(INT(1), *N)) && (!lQuery)) {
        *INFO = -8;
    }
}


void
ORGQR(INT *M, INT *N, INT *K, FLOAT *A, INT *LDA, FLOAT *TAU,
      FLOAT *WORK, INT *LWORK, INT *INFO)
{
    orgqrErrorCheck(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(ORGQR_NAME, INFO, strlen(ORGQR_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifdef USE_ORGQR_REF
    // call LAPACK implementation

    ORGQR_REF(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);

#   elif defined(CMP_W_ORGQR_REF)
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

    lapack::orgqr(*K, __A, __tau, __work);
    ORGQR_REF(M, N, K, A, LDA, TAU, WORK, LWORK, INFO);

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

    lapack::orgqr(*K, _A, _tau, _work);
#   endif

}

} // extern "C"

