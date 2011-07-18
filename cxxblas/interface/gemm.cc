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

extern "C" {

    void xerbla_(const char* srname, int* info);

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_TRANSPOSE    {CblasNoTrans=111, CblasTrans=112,
                             CblasConjTrans=113, CblasConjNoTrans=114};

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CBLAS_ALPHA;
    typedef CBLAS_FLOAT                 CBLAS_BETA;
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef const CBLAS_FLOAT *         CBLAS_ALPHA;
    typedef const CBLAS_FLOAT *         CBLAS_BETA;
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_TRANSPOSE _transA, enum CBLAS_TRANSPOSE _transB,
               CBLAS_INT m, CBLAS_INT n, CBLAS_INT k,
               CBLAS_ALPHA _alpha,
               const CBLAS_FLOAT  *_A, CBLAS_INT ldA,
               const CBLAS_FLOAT  *_B, CBLAS_INT ldB,
               CBLAS_BETA _beta,
               CBLAS_FLOAT  *_C, CBLAS_INT ldC)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;
        Transpose transA = NoTrans;
        CBLAS_INT numRowsA = m;
        CBLAS_INT numColsA = k;
        if (_transA==CblasTrans) {
            transA = Trans;
            numRowsA = k;
            numColsA = m;
        }
        if (_transA==CblasConjTrans) {
            transA = ConjTrans;
            numRowsA = k;
            numColsA = m;
        }
        if (_transA==CblasConjNoTrans) {
            transA = Conj;
            numRowsA = m;
            numColsA = k;
        }

        Transpose transB = NoTrans;
        CBLAS_INT numRowsB = k;
        CBLAS_INT numColsB = n;
        if (_transB==CblasTrans) {
            transB = Trans;
            numRowsB = n;
            numColsB = k;
        }
        if (_transB==CblasConjTrans) {
            transB = ConjTrans;
            numRowsB = n;
            numColsB = k;
        }
        if (_transB==CblasConjNoTrans) {
            transB = Conj;
            numRowsB = k;
            numColsB = n;
        }

        const CXXBLAS_FLOAT *A = reinterpret_cast<const CXXBLAS_FLOAT *>(_A);
        const CXXBLAS_FLOAT *B = reinterpret_cast<const CXXBLAS_FLOAT *>(_B);
        CXXBLAS_FLOAT *C = reinterpret_cast<CXXBLAS_FLOAT *>(_C);

#   ifndef COMPLEX_FLOAT1
        CXXBLAS_FLOAT alpha = _alpha;
        CXXBLAS_FLOAT beta = _beta;
#   else
        CXXBLAS_FLOAT alpha(_alpha[0], _alpha[1]);
        CXXBLAS_FLOAT beta(_beta[0], _beta[1]);
#   endif

        CXXBLAS_FLOAT *__C = C;
        CBLAS_INT      __ldC = ldC;
        StorageOrder   __order = order;

#   ifdef TEST_ROW_MAJOR
        switchFullStorageOrder(order, numRowsA, numColsA, A, ldA);
        switchFullStorageOrder(order, numRowsB, numColsB, B, ldB);

        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocateFullStorage(__order, m, n, __C, __ldC);
        switchFullStorageOrder(order, m, n, C, ldC, __C, __ldC);
#   endif

        cxxblas::gemm(__order, transA, transB,
                      m, n, k,
                      alpha,
                      A, ldA,
                      B, ldB,
                      beta,
                      __C, __ldC);

#   ifdef TEST_ROW_MAJOR
        releaseStorage(A);
        releaseStorage(B);
        switchFullStorageOrder(__order, m, n, __C, __ldC, C, ldC);
        releaseStorage(__C);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const char *_transA, const char *_transB,
              const CBLAS_INT *_m, const CBLAS_INT *_n, const CBLAS_INT *_k,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *A, const CBLAS_INT *_ldA,
              const CBLAS_FLOAT *B, const CBLAS_INT *_ldB,
              const CBLAS_FLOAT *_beta,
              CBLAS_FLOAT *C, const CBLAS_INT *_ldC)
    {
        CBLAS_TRANSPOSE transA, transB;
        CBLAS_INT numRowsA = -1;
        CBLAS_INT numRowsB = -1;

        bool checkTransA = false;
        bool checkTransB = false;

        CBLAS_INT m       = *_m;
        CBLAS_INT n       = *_n;
        CBLAS_INT k       = *_k;
        CBLAS_INT ldA     = *_ldA;
        CBLAS_INT ldB     = *_ldB;
        CBLAS_INT ldC     = *_ldC;

        if ((*_transA=='N') || (*_transA=='n')) {
            transA = CblasNoTrans;
            checkTransA = true;
            numRowsA = m;
        }
        if ((*_transA=='T') || (*_transA=='t')) {
            transA = CblasTrans;
            checkTransA = true;
            numRowsA = k;
        }
        if ((*_transA=='C') || (*_transA=='c')) {
            transA = CblasConjTrans;
            checkTransA = true;
            numRowsA = k;
        }
        if ((*_transA=='R') || (*_transA=='r')) {
            transA = CblasConjNoTrans;
            checkTransA = true;
            numRowsA = m;
        }

        if ((*_transB=='N') || (*_transB=='n')) {
            transB = CblasNoTrans;
            checkTransB = true;
            numRowsB = k;
        }
        if ((*_transB=='T') || (*_transB=='t')) {
            transB = CblasTrans;
            checkTransB = true;
            numRowsB = n;
        }
        if ((*_transB=='C') || (*_transB=='c')) {
            transB = CblasConjTrans;
            checkTransB = true;
            numRowsB = n;
        }
        if ((*_transB=='R') || (*_transB=='r')) {
            transB = CblasConjNoTrans;
            checkTransB = true;
            numRowsB = k;
        }

#   ifndef COMPLEX_FLOAT1
        CBLAS_ALPHA alpha  = *_alpha;
        CBLAS_BETA beta  = *_beta;
#   else
        CBLAS_ALPHA alpha = _alpha;
        CBLAS_BETA beta = _beta;
#   endif

        CBLAS_INT info = 0;
        if (ldC<std::max(1, m)) {
            info = 13;
        }
        if (ldB<std::max(1, numRowsB)) {
            info = 10;
        }
        if (ldA<std::max(1, numRowsA)) {
            info = 8;
        }
        if (k<0) {
            info = 5;
        }
        if (n<0) {
            info = 4;
        }
        if (m<0) {
            info = 3;
        }
        if (!checkTransB) {
            info = 2;
        }
        if (!checkTransA) {
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
        CBLAS_NAME(CblasColMajor, transA, transB,
                   m, n, k,
                   alpha,
                   A, ldA,
                   B, ldB,
                   beta,
                   C, ldC);
    }
#endif // CREATE_BLAS

} // extern "C"