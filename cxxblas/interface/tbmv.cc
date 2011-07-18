#include <cxxblas/cxxblas.cxx>
#include <cxxblas/interface/aux.h>

using cxxblas::StorageOrder;
using cxxblas::ColMajor;
using cxxblas::RowMajor;

using cxxblas::Transpose;
using cxxblas::NoTrans;
using cxxblas::Trans;
using cxxblas::Conj;
using cxxblas::ConjTrans;

using cxxblas::StorageUpLo;
using cxxblas::Upper;
using cxxblas::Lower;

using cxxblas::Diag;
using cxxblas::Unit;
using cxxblas::NonUnit;


extern "C" {

    void xerbla_(const char* srname, int* info);

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_UPLO         {CblasUpper=121, CblasLower=122};
    enum CBLAS_TRANSPOSE    {CblasNoTrans=111, CblasTrans=112,
                             CblasConjTrans=113, CblasConjNoTrans=114};
    enum CBLAS_DIAG         {CblasNonUnit=131, CblasUnit=132};


    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_UPLO _upLo,
               enum CBLAS_TRANSPOSE _transA,
               enum CBLAS_DIAG _diag,
               CBLAS_INT n, CBLAS_INT k,
               const CBLAS_FLOAT  *_A, CBLAS_INT ldA,
               CBLAS_FLOAT  *_x, CBLAS_INT incX)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;
        StorageUpLo upLo = Lower;
        if (_upLo==CblasUpper) {
            upLo = Upper;
        }

        Transpose transA = NoTrans;
        if (_transA==CblasTrans) {
            transA = Trans;
        }
        if (_transA==CblasConjTrans) {
            transA = ConjTrans;
        }
        if (_transA==CblasConjNoTrans) {
            transA = Conj;
        }

        Diag diag = NonUnit;
        if (_diag==CblasUnit) {
            diag = Unit;
        }

        const CXXBLAS_FLOAT *A = reinterpret_cast<const CXXBLAS_FLOAT *>(_A);
        CXXBLAS_FLOAT *x = reinterpret_cast<CXXBLAS_FLOAT *>(_x);

#   ifdef TEST_ROW_MAJOR
        switchBandStorageOrder(order, upLo, n, k, A, ldA);
        order = (order==ColMajor) ? RowMajor : ColMajor;
#   endif

        cxxblas::tbmv(order, upLo, transA, diag,
                      n, k,
                      A, ldA,
                      x, incX);

#   ifdef TEST_ROW_MAJOR
        releaseStorage(A);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const char *_upLo, const char *_transA, const char *_diag,
              const CBLAS_INT *_n, const CBLAS_INT *_k,
              const CBLAS_FLOAT *A, const CBLAS_INT *_ldA,
              CBLAS_FLOAT *x, const CBLAS_INT *_incX)
    {
        bool checkUpLo = false,
             checkDiag = false,
             checkTrans = false;
        
        CBLAS_UPLO upLo = CblasUpper;
        if ((*_upLo=='L') || (*_upLo=='l')) {
            upLo = CblasLower;
            checkUpLo = true;
        }
        if ((*_upLo=='U') || (*_upLo=='u')) {
            checkUpLo = true;
        }

        CBLAS_TRANSPOSE transA = CblasConjNoTrans;
        if ((*_transA=='N') || (*_transA=='n')) {
            transA = CblasNoTrans;
            checkTrans = true;
        }
        if ((*_transA=='T') || (*_transA=='t')) {
            transA = CblasTrans;
            checkTrans = true;
        }
        if ((*_transA=='C') || (*_transA=='c')) {
            transA = CblasConjTrans;
            checkTrans = true;
        }
        if ((*_transA=='R') || (*_transA=='r')) {
            transA = CblasConjNoTrans;
            checkTrans = true;
        }

        CBLAS_DIAG diag = CblasUnit;
        if ((*_diag=='N') || (*_diag=='n')) {
            diag = CblasNonUnit;
            checkDiag = true;
        }
        if ((*_diag=='U') || (*_diag=='u')) {
            checkDiag = true;
        }

        CBLAS_INT n       = *_n;
        CBLAS_INT k       = *_k;
        CBLAS_INT ldA     = *_ldA;
        CBLAS_INT incX    = *_incX;

        CBLAS_INT info = 0;
        if (incX==0) {
            info = 9; 
        }
        if (ldA<k+1) {
            info = 7;
        }
        if (k<0) {
            info = 5;
        }
        if (n<0) {
            info = 4;
        }
        if (!checkDiag) {
            info = 3;
        }
        if (!checkTrans) {
            info = 2;
        }
        if (!checkUpLo) {
            info = 1;
        }
        if (info!=0) {
            char blasName[6];
            strncpy(blasName, BLAS_NAME_STR, 6);
            for (int i=0; i<6; ++i) {
                blasName[i] = std::toupper(blasName[i]);
            }
            xerbla_(blasName, &info);
          return;
        }

        // the blas interface calls the cblas interface
        // so any blas-test will also test the cblas-interface
        CBLAS_NAME(CblasColMajor, upLo, transA, diag,
                   n, k,
                   A, ldA,
                   x, incX);
    }
#endif // CREATE_BLAS

} // extern "C"