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
    typedef CBLAS_FLOAT                 CBLAS_BETA;
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef const CBLAS_FLOAT *         CBLAS_ALPHA;
    typedef const CBLAS_FLOAT *         CBLAS_BETA;
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_UPLO _upLoC, enum CBLAS_TRANSPOSE _transA,
               CBLAS_INT n, CBLAS_INT k,
               CBLAS_ALPHA _alpha,
               const CBLAS_FLOAT  *_A, CBLAS_INT ldA,
               CBLAS_BETA _beta,
               CBLAS_FLOAT  *_C, CBLAS_INT ldC)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;

        StorageUpLo upLoC = (_upLoC==CblasUpper) ? Upper : Lower;

        Transpose transA = NoTrans;
        CBLAS_INT numRowsA = n;
        CBLAS_INT numColsA = k;

        if (_transA==CblasTrans) {
            transA = Trans;
            numRowsA = k;
            numColsA = n;
        }
        if (_transA==CblasConjTrans) {
            transA = ConjTrans;
            numRowsA = k;
            numColsA = n;
        }
        if (_transA==CblasConjNoTrans) {
            transA = Conj;
            numRowsA = n;
            numColsA = k;
        }

        const CXXBLAS_FLOAT *A = reinterpret_cast<const CXXBLAS_FLOAT *>(_A);
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

        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocateFullStorage(__order, n, n, __C, __ldC);
        switchFullStorageOrder(order, n, n, C, ldC, __C, __ldC);
#   endif

        cxxblas::syrk(__order, upLoC, transA,
                      n, k,
                      alpha,
                      A, ldA,
                      beta,
                      __C, __ldC);

#   ifdef TEST_ROW_MAJOR
        releaseStorage(A);

        switchFullStorageOrder(__order, n, n, __C, __ldC, C, ldC);
        releaseStorage(__C);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS


#ifdef CREATE_BLAS
    void
    BLAS_NAME(const char *_upLoC, const char *_transA,
              const CBLAS_INT *_n, const CBLAS_INT *_k,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *A, const CBLAS_INT *_ldA,
              const CBLAS_FLOAT *_beta,
              CBLAS_FLOAT *C, const CBLAS_INT *_ldC)
    {
        CBLAS_UPLO          upLoC;
        CBLAS_TRANSPOSE     transA;
        CBLAS_INT           numRowsA  = -1;

        bool checkUpLoC = false;
        bool checkTransA = false;

        CBLAS_INT n       = *_n;
        CBLAS_INT k       = *_k;
        CBLAS_INT ldA     = *_ldA;
        CBLAS_INT ldC     = *_ldC;

        if ((*_upLoC=='L') || (*_upLoC=='l')) {
            upLoC = CblasLower;
            checkUpLoC = true;
        }
        if ((*_upLoC=='U') || (*_upLoC=='u')) {
            upLoC = CblasUpper;
            checkUpLoC = true;
        }

        if ((*_transA=='N') || (*_transA=='n')) {
            transA = CblasNoTrans;
            numRowsA = n;
            checkTransA = true;
        }
        if ((*_transA=='T') || (*_transA=='t')) {
            transA = CblasTrans;
            numRowsA = k;
            checkTransA = true;
        }
        if ((*_transA=='C') || (*_transA=='c')) {
            transA = CblasConjTrans;
            numRowsA = k;
#       ifndef COMPLEX_FLOAT1
            checkTransA = true;
#       else
            checkTransA = false;
#       endif
        }
        if ((*_transA=='R') || (*_transA=='r')) {
            transA = CblasConjNoTrans;
            numRowsA = n;
#       ifndef COMPLEX_FLOAT1
            checkTransA = true;
#       else
            checkTransA = false;
#       endif
        }


#   ifndef COMPLEX_FLOAT1
        CBLAS_ALPHA alpha  = *_alpha;
        CBLAS_BETA beta  = *_beta;
#   else
        CBLAS_ALPHA alpha = _alpha;
        CBLAS_ALPHA beta = _beta;
#   endif

        CBLAS_INT info = 0;
        if (ldC<std::max(1, n)) {
            info = 10;
        }
        if (ldA<std::max(1, numRowsA)) {
            info = 7;
        }
        if (k<0) {
            info = 4;
        }
        if (n<0) {
            info = 3;
        }
        if (!checkTransA) {
            info = 2;
        }
        if (!checkUpLoC) {
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
        CBLAS_NAME(CblasColMajor, upLoC, transA,
                   n, k,
                   alpha,
                   A, ldA,
                   beta,
                   C, ldC);
    }
#endif // CREATE_BLAS

} // extern "C"