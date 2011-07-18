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
               enum CBLAS_TRANSPOSE _transA,
               CBLAS_INT m, CBLAS_INT n,
               CBLAS_ALPHA _alpha,
               const CBLAS_FLOAT  *_A, CBLAS_INT ldA,
               const CBLAS_FLOAT  *_x, CBLAS_INT incX,
               CBLAS_BETA _beta,
               CBLAS_FLOAT  *_y, CBLAS_INT incY)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;
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
        const CXXBLAS_FLOAT *x = reinterpret_cast<const CXXBLAS_FLOAT *>(_x);
        CXXBLAS_FLOAT *y = reinterpret_cast<CXXBLAS_FLOAT *>(_y);

#   ifndef COMPLEX_FLOAT1
        CXXBLAS_FLOAT alpha = _alpha;
        CXXBLAS_FLOAT beta = _beta;
#   else
        CXXBLAS_FLOAT alpha(_alpha[0], _alpha[1]);
        CXXBLAS_FLOAT beta(_beta[0], _beta[1]);
#   endif

#   ifdef TEST_ROW_MAJOR
        switchFullStorageOrder(order, m, n, A, ldA);
        order = (order==ColMajor) ? RowMajor : ColMajor;
#   endif

        cxxblas::gemv(order, transA,
                      m, n,
                      alpha,
                      A, ldA,
                      x, incX,
                      beta,
                      y, incY);

#   ifdef TEST_ROW_MAJOR
        releaseStorage(A);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const char *_transA,
              const CBLAS_INT *_m, const CBLAS_INT *_n,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *A, const CBLAS_INT *_ldA,
              const CBLAS_FLOAT *x, const CBLAS_INT *_incX,
              const CBLAS_FLOAT *_beta,
              CBLAS_FLOAT *y, const CBLAS_INT *_incY)
    {
        CBLAS_TRANSPOSE transA;

        bool checkTransA = false;
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

        CBLAS_INT m       = *_m;
        CBLAS_INT n       = *_n;
        CBLAS_INT ldA     = *_ldA;
        CBLAS_INT incX    = *_incX;
        CBLAS_INT incY    = *_incY;

#   ifndef COMPLEX_FLOAT1
        CBLAS_ALPHA alpha  = *_alpha;
        CBLAS_BETA beta  = *_beta;
#   else
        CBLAS_ALPHA alpha = _alpha;
        CBLAS_BETA beta = _beta;
#   endif

        CBLAS_INT info = 0;
        if (incY==0) {
            info = 11;
        }
        if (incX==0) {
            info = 8; 
        }
        if (ldA<std::max(1, m)) {
            info = 6;
        }
        if (n<0) {
            info = 3;
        }
        if (m<0) {
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
        CBLAS_NAME(CblasColMajor, transA,
                   m, n,
                   alpha,
                   A, ldA,
                   x, incX,
                   beta,
                   y, incY);
    }
#endif // CREATE_BLAS

} // extern "C"