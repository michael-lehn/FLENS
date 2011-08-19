#include <flens/lapack/interface/interface.h>

#ifdef SINGLE
#   define GESV         sgesv_
#   define GESV_REF     sgesv
#   define GESV_NAME    "SGESV"
#elif DOUBLE
#   define GESV         dgesv_
#   define GESV_REF     dgesv
#   define GESV_NAME    "DGESV"
#elif COMPLEX_SINGLE
#   define GESV         cgesv_
#   define GESV_REF     cgesv
#   define GESV_NAME    "CGESV"
#elif COMPLEX_DOUBLE
#   define GESV         zgesv_
#   define GESV_REF     zgesv
#   define GESV_NAME    "ZGESV"
#endif

using namespace flens;

extern "C" {


void
GESV_REF(INT *N, INT *NRHS,
         FLOAT *A, INT *LDA,
         INT *IPIV,
         FLOAT *B, INT *LDB,
         INT *INFO);


void
gesvErrorCheck(INT *N, INT *NRHS,
               FLOAT *A, INT *LDA,
               INT *IPIV,
               FLOAT *B, INT *LDB,
               INT *INFO)
{
    if (*N<0) {
        *INFO = -1;
        return;
    }
    if (*NRHS<0) {
        *INFO = -2;
        return;
    }
    if (*LDA<std::max(INT(1), *N)) {
        *INFO = -4;
        return;
    }
    if (*LDB<std::max(INT(1), *N)) {
        *INFO = -7;
        return;
    }
}


void
GESV(INT *N, INT *NRHS,
     FLOAT *A, INT *LDA,
     INT *IPIV,
     FLOAT *B, INT *LDB,
     INT *INFO)
{
    gesvErrorCheck(N, NRHS, A, LDA, IPIV, B, LDB, INFO);

    if (*INFO<0) {
        *INFO = -(*INFO);
        XERBLA(GESV_NAME, INFO, strlen(GESV_NAME));
        *INFO = -(*INFO);
        return;
    }


#   ifdef USE_GESV_REF
    // call LAPACK implementation

    GESV_REF(N, NRHS, A, LDA, IPIV, B, LDB, INFO);

#   elif defined(CMP_W_GESV_REF)
    // call LAPACK and FLENS implementation, compare results

    if ((*N==0) || (*NRHS==0)) {
        *INFO = 0;
        return;
    }

    GeMatrix<FSV>        _A = FSV(*N, *N, A, *LDA);
    DenseVector<IAV>     _piv = IAV(*N, IPIV, INT(1));
    GeMatrix<FSV>        _B = FSV(*N, *NRHS, B, *LDB);

    GeMatrix<FS>         __A = _A;
    DenseVector<IA>      __piv = _piv;
    GeMatrix<FS>         __B = _B;


    *INFO = lapack::sv(__A, __piv, __B);
    GESV_REF(N, NRHS, A, LDA, IPIV, B, LDB, INFO);

    if (isDifferent(_A,__A)) {
        std::cerr << "_A = " << _A << std::endl;
        std::cerr << "__A = " << __A << std::endl;
        assert(0);
    }

    if (isDifferent(_B,__B)) {
        std::cerr << "_B = " << _B << std::endl;
        std::cerr << "__B = " << __B << std::endl;
        assert(0);
    }

    if (isDifferent(_piv,__piv)) {
        std::cerr << "_piv = " << _A << std::endl;
        std::cerr << "__piv = " << __A << std::endl;
        assert(0);
    }

#   else
    // call FLENS implementation

    GeMatrix<FSV>        _A = FSV(*N, *N, A, *LDA);
    GeMatrix<FSV>        _B = FSV(*N, *NRHS, B, *LDB);
    DenseVector<IAV>     _piv = IAV(*N, IPIV, INT(1));

    *INFO = lapack::sv(_A, _piv, _B);
#   endif
}

} // extern "C"

