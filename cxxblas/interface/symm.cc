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
               enum CBLAS_SIDE _sideA,  enum CBLAS_UPLO _upLoA,
               CBLAS_INT m, CBLAS_INT n,
               CBLAS_ALPHA _alpha,
               const CBLAS_FLOAT  *_A, CBLAS_INT ldA,
               const CBLAS_FLOAT  *_B, CBLAS_INT ldB,
               CBLAS_BETA _beta,
               CBLAS_FLOAT  *_C, CBLAS_INT ldC)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;

        Side        sideA = (_sideA==CblasLeft)  ? Left : Right;
        StorageUpLo upLoA = (_upLoA==CblasUpper) ? Upper : Lower;

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
        CBLAS_INT   dimA =  (sideA==Left)        ? m : n;

        switchFullStorageOrder(order, dimA, dimA, A, ldA);
        switchFullStorageOrder(order, m, n, B, ldB);

        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocateFullStorage(__order, m, n, __C, __ldC);
        switchFullStorageOrder(order, m, n, C, ldC, __C, __ldC);
#   endif

        cxxblas::symm(__order, sideA, upLoA,
                      m, n,
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
    BLAS_NAME(const char *_sideA, const char *_upLoA,
              const CBLAS_INT *_m, const CBLAS_INT *_n,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *A, const CBLAS_INT *_ldA,
              const CBLAS_FLOAT *B, const CBLAS_INT *_ldB,
              const CBLAS_FLOAT *_beta,
              CBLAS_FLOAT *C, const CBLAS_INT *_ldC)
    {
        CBLAS_SIDE          sideA;
        CBLAS_UPLO          upLoA;
        CBLAS_INT           dimA  = -1;

        bool checkSideA = false;
        bool checkUpLoA = false;

        CBLAS_INT m       = *_m;
        CBLAS_INT n       = *_n;
        CBLAS_INT ldA     = *_ldA;
        CBLAS_INT ldB     = *_ldB;
        CBLAS_INT ldC     = *_ldC;

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

#   ifndef COMPLEX_FLOAT1
        CBLAS_ALPHA alpha  = *_alpha;
        CBLAS_BETA beta  = *_beta;
#   else
        CBLAS_ALPHA alpha = _alpha;
        CBLAS_ALPHA beta = _beta;
#   endif

        CBLAS_INT info = 0;
        if (ldC<std::max(1, m)) {
            info = 12;
        }
        if (ldB<std::max(1, m)) {
            info = 9;
        }
        if (ldA<std::max(1, dimA)) {
            info = 7;
        }
        if (n<0) {
            info = 4;
        }
        if (m<0) {
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
        CBLAS_NAME(CblasColMajor, sideA, upLoA,
                   m, n,
                   alpha,
                   A, ldA,
                   B, ldB,
                   beta,
                   C, ldC);
    }
#endif // CREATE_BLAS

} // extern "C"