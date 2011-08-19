#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GETRF         sgetrf_
#   define GETRF_REF     sgetrf
#   define GETRF_NAME    "SGETRF"
#elif DOUBLE
#   define GETRF         dgetrf_
#   define GETRF_REF     dgetrf
#   define GETRF_NAME    "DGETRF"
#elif COMPLEX_SINGLE
#   define GETRF         cgetrf_
#   define GETRF_REF     cgetrf
#   define GETRF_NAME    "CGETRF"
#elif COMPLEX_DOUBLE
#   define GETRF         zgetrf_
#   define GETRF_REF     zgetrf
#   define GETRF_NAME    "ZGETRF"
#endif


using namespace flens;

extern "C" {

void
GETRF_REF(INT *M, INT *N, FLOAT *A, INT *LDA, INT *IPIV, INT *INFO);


void
getrfErrorCheck(INT *M, INT *N, FLOAT *A, INT *LDA, INT *IPIV, INT *INFO)
{
    if (*M<0) {
        *INFO = -1;
    } else if (*N<0) {
        *INFO = -2;
    } else if (*LDA<std::max(INT(1), *M)) {
        *INFO = -4;
        return;
    }
}


void
GETRF(INT *M, INT *N, FLOAT *A, INT *LDA, INT *IPIV, INT *INFO)
{
    getrfErrorCheck(M, N, A, LDA, IPIV, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GETRF_NAME, INFO, strlen(GETRF_NAME));
        *INFO = -(*INFO);
        return;
    }


#   ifdef USE_GETRF_REF
    // call LAPACK implementation

    GETRF_REF(M, N, A, LDA, IPIV, INFO);

#   elif defined(CMP_W_GETRF_REF)
    // call LAPACK and FLENS implementation, compare results

    if ((*M==0) || (*N==0)) {
        *INFO = 0;
        return;
    }

    GeMatrix<FSV>        _A = FSV(*M, *N, A, *LDA);
    DenseVector<IAV>     _piv = IAV(*M, IPIV, INT(1));

    GeMatrix<FS>         __A = _A;
    DenseVector<IA>      __piv = _piv;


    INT info = lapack::trf(__A, __piv);
    GETRF_REF(M, N, A, LDA, IPIV, INFO);

    if (info!=*INFO) {
        std::cerr << "info = " << info << ", *INFO = " << *INFO << std::endl;
    }

    if (isDifferent(_A,__A) || isDifferent(_piv,__piv)) {
        assert(0);
    }

#   else
    // call FLENS implementation

    GeMatrix<FSV>        _A = FSV(*M, *N, A, *LDA);
    DenseVector<IAV>     _piv = IAV(*M, IPIV, INT(1));

    *INFO = lapack::trf(_A, _piv);
#   endif
}

} // extern "C"

