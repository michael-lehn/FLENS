#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GEQRF         sgeqrf_
#   define GEQRF_REF     sgeqrf
#   define GEQRF_NAME    "SGEQRF"
#elif DOUBLE
#   define GEQRF         dgeqrf_
#   define GEQRF_REF     dgeqrf
#   define GEQRF_NAME    "DGEQRF"
#elif COMPLEX_SINGLE
#   define GEQRF         cgeqrf_
#   define GEQRF_REF     cgeqrf
#   define GEQRF_NAME    "CGEQRF"
#elif COMPLEX_DOUBLE
#   define GEQRF         zgeqrf_
#   define GEQRF_REF     zgeqrf
#   define GEQRF_NAME    "ZGEQRF"
#endif

using namespace flens;

extern "C" {

void
GEQRF_REF(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
          FLOAT *WORK, INT *LWORK, INT *INFO);

void
geqrfErrorCheck(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
                FLOAT *WORK, INT *LWORK, INT *INFO)
{
    bool lQuery = (*LWORK==-1);
    
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -4;
    } else if ((*LWORK<std::max(INT(1), *N)) && (!lQuery)) {
        *INFO = -7;
    }
}


void
GEQRF(INT *M, INT *N, FLOAT *A, INT *LDA, FLOAT *TAU,
      FLOAT *WORK, INT *LWORK, INT *INFO)
{
#   ifdef DEBUG_INTERFACE
    std::cerr << GEQRF_NAME << std::endl;
#   endif

    geqrfErrorCheck(M, N, A, LDA, TAU, WORK, LWORK, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GEQRF_NAME, INFO, strlen(GEQRF_NAME));
        *INFO = -(*INFO);
        return;
    }


#   ifdef USE_GEQRF_REF
    // call LAPACK implementation

    GEQRF_REF(M, N, A, LDA, TAU, WORK, LWORK, INFO);

#   elif defined(CMP_W_GEQRF_REF)
    // call LAPACK and FLENS implementation, compare results

    if ((*M==0) || (*N==0)) {
        *INFO = 0;
        return;
    }

    INT k = std::min(*M, *N);

    GeMatrix<FSV>        _A = FSV(*M, *N, A, *LDA);
    DenseVector<AV>      _tau = AV(k, TAU, INT(1));
    DenseVector<AV>      _work = AV(*LWORK, WORK, INT(1));

    GeMatrix<FS>         __A = _A;
    DenseVector<AR>      __tau = _tau;
    DenseVector<AR>      __work = _work;


    lapack::qrf(__A, __tau, __work);
    GEQRF_REF(M, N, A, LDA, TAU, WORK, LWORK, INFO);

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
        //std::cerr << "_work = " << _work << std::endl;
        //std::cerr << "__work = " << __work << std::endl;
        //assert(0);
    }

#   else
    // call FLENS implementation

    GeMatrix<FSV>        _A = FSV(*M, *N, A, *LDA);
    DenseVector<AV>      _tau = AV(k, TAU, INT(1));
    DenseVector<AV>      _work = AV(*LWORK, WORK, INT(1));

    lapack::qrf(__A, __tau, __work);
#   endif


}

} // extern "C"

