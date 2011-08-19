#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GETF2         sgetf2_
#   define GETF2_REF     sgetf2
#   define GETF2_NAME    "SGETF2"
#elif DOUBLE
#   define GETF2         dgetf2_
#   define GETF2_REF     dgetf2
#   define GETF2_NAME    "DGETF2"
#elif COMPLEX_SINGLE
#   define GETF2         cgetf2_
#   define GETF2_REF     cgetf2
#   define GETF2_NAME    "CGETF2"
#elif COMPLEX_DOUBLE
#   define GETF2         zgetf2_
#   define GETF2_REF     zgetf2
#   define GETF2_NAME    "ZGETF2"
#endif

using namespace flens;

extern "C" {

void
GETF2_REF(INT *M, INT *N, FLOAT *A, INT *LDA, INT *IPIV, INT *INFO);


void
getf2ErrorCheck(INT *M, INT *N, FLOAT *A, INT *LDA, INT *IPIV, INT *INFO)
{
    if (*M<0) {
        *INFO = -1;
        return;
    }
    if (*N<0) {
        *INFO = -2;
        return;
    }
    if (*LDA<std::max(INT(1), *M)) {
        *INFO = -4;
        return;
    }
}


void
GETF2(INT *M, INT *N, FLOAT *A, INT *LDA, INT *IPIV, INT *INFO)
{
    getf2ErrorCheck(M, N, A, LDA, IPIV, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GETF2_NAME, INFO, strlen(GETF2_NAME));
        *INFO = -(*INFO);
        return;
    }

#   ifdef USE_GETRF_REF
    // call LAPACK implementation

    GETF2_REF(M, N, A, LDA, IPIV, INFO);

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


    INT info = lapack::tf2(__A, __piv);
    GETF2_REF(M, N, A, LDA, IPIV, INFO);

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

    *INFO = lapack::tf2(_A, _piv);
#   endif
}

} // extern "C"

