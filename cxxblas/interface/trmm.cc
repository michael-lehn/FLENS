#include <cxxblas/cxxblas.cxx>
#include <cxxblas/interface/aux.h>
#include <complex>

using cxxblas::StorageOrder;
using cxxblas::ColMajor;
using cxxblas::RowMajor;
using cxxblas::Transpose;
using cxxblas::NoTrans;
using cxxblas::Trans;
using cxxblas::Conj;
using cxxblas::ConjTrans;
using cxxblas::Side;
using cxxblas::Left;
using cxxblas::Right;
using cxxblas::Diag;
using cxxblas::Unit;
using cxxblas::NonUnit;


extern "C" {

    void xerbla_(const char* srname, int* info);

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_TRANSPOSE    {CblasNoTrans=111, CblasTrans=112,
                             CblasConjTrans=113, CblasConjNoTrans=114};
    enum CBLAS_SIDE         {CblasLeft=141, CblasRight=142};
    enum CBLAS_UPLO         {CblasUpper=121, CblasLower=122};
    enum CBLAS_DIAG         {CblasNonUnit=131, CblasUnit=132};

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CBLAS_ALPHA;
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef const CBLAS_FLOAT *         CBLAS_ALPHA;
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_SIDE _sideA,  enum CBLAS_UPLO _upLoA,
               enum CBLAS_TRANSPOSE _transA, enum CBLAS_DIAG _diagA,
               CBLAS_INT m, CBLAS_INT n,
               CBLAS_ALPHA _alpha,
               const CBLAS_FLOAT  *_A, CBLAS_INT ldA,
               CBLAS_FLOAT  *_B, CBLAS_INT ldB)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;

        Side        sideA = (_sideA==CblasLeft)  ? Left : Right;
        StorageUpLo upLoA = (_upLoA==CblasUpper) ? Upper : Lower;
        Diag        diagA = (_diagA==CblasUnit)  ? Unit : NonUnit;

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

        const CXXBLAS_FLOAT *A = reinterpret_cast<const CXXBLAS_FLOAT *>(_A);
        CXXBLAS_FLOAT *B = reinterpret_cast<CXXBLAS_FLOAT *>(_B);

#   ifndef COMPLEX_FLOAT1
        CXXBLAS_FLOAT alpha = _alpha;
#   else
        CXXBLAS_FLOAT alpha(_alpha[0], _alpha[1]);
#   endif

        CXXBLAS_FLOAT *__B = B;
        CBLAS_INT      __ldB = ldB;
        StorageOrder   __order = order;

#   ifdef TEST_ROW_MAJOR
        CBLAS_INT   dimA =  (sideA==Left)        ? m : n;

        switchFullStorageOrder(order, dimA, dimA, A, ldA);

        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocateFullStorage(__order, m, n, __B, __ldB);
        switchFullStorageOrder(order, m, n, B, ldB, __B, __ldB);
#   endif

        cxxblas::trmm(__order, sideA, upLoA, transA, diagA,
                      m, n,
                      alpha,
                      A, ldA,
                      __B, __ldB);

#   ifdef TEST_ROW_MAJOR
        releaseStorage(A);
        switchFullStorageOrder(__order, m, n, __B, __ldB, B, ldB);
        releaseStorage(__B);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const char *_sideA, const char *_upLoA,
              const char *_transA, const char *_diagA,
              const CBLAS_INT *_m, const CBLAS_INT *_n,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *A, const CBLAS_INT *_ldA,
              CBLAS_FLOAT *B, const CBLAS_INT *_ldB)
    {
        CBLAS_SIDE          sideA;
        CBLAS_UPLO          upLoA;
        CBLAS_TRANSPOSE     transA;
        CBLAS_DIAG          diagA;
        CBLAS_INT           dimA  = -1;

        bool checkSideA = false;
        bool checkUpLoA = false;
        bool checkTransA = false;
        bool checkDiagA = false;

        CBLAS_INT m       = *_m;
        CBLAS_INT n       = *_n;
        CBLAS_INT ldA     = *_ldA;
        CBLAS_INT ldB     = *_ldB;

        if (*_sideA=='L') {
            sideA = CblasLeft;
            dimA = m;
            checkSideA = true;
        }
        if (*_sideA=='R') {
            sideA = CblasRight;
            dimA = n;
            checkSideA = true;
        }

        if ((*_upLoA=='L') || (*_upLoA=='l')) {
            upLoA = CblasLower;
            checkUpLoA = true;
        }
        if ((*_upLoA=='U') || (*_upLoA=='u')) {
            upLoA = CblasUpper;
            checkUpLoA = true;
        }

        if ((*_transA=='N') || (*_transA=='n')) {
            transA = CblasNoTrans;
            checkTransA = true;
        }
        if ((*_transA=='T') || (*_transA=='t')) {
            transA = CblasTrans;
            checkTransA = true;
        }
        if ((*_transA=='C') || (*_transA=='c')) {
            transA = CblasConjTrans;
            checkTransA = true;
        }
        if ((*_transA=='R') || (*_transA=='r')) {
            transA = CblasConjNoTrans;
            checkTransA = true;
        }

        if ((*_diagA=='N') || (*_diagA=='n')) {
            diagA = CblasNonUnit;
            checkDiagA = true;
        }
        if ((*_diagA=='U') || (*_diagA=='u')) {
            diagA = CblasUnit;
            checkDiagA = true;
        }


#   ifndef COMPLEX_FLOAT1
        CBLAS_ALPHA alpha  = *_alpha;
#   else
        CBLAS_ALPHA alpha = _alpha;
#   endif

        CBLAS_INT info = 0;
        if (ldB<std::max(1, m)) {
            info = 11;
        }
        if (ldA<std::max(1, dimA)) {
            info = 9;
        }
        if (n<0) {
            info = 6;
        }
        if (m<0) {
            info = 5;
        }
        if (!checkDiagA) {
            info = 4;
        }
        if (!checkTransA) {
            info = 3;
        }
        if (!checkUpLoA) {
            info = 2;
        }
        if (!checkSideA) {
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
        CBLAS_NAME(CblasColMajor, sideA, upLoA, transA, diagA,
                   m, n,
                   alpha,
                   A, ldA,
                   B, ldB);
    }
#endif // CREATE_BLAS

} // extern "C"