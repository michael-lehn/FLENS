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

extern "C" {

    void xerbla_(const char* srname, int* info);

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_TRANSPOSE    {CblasNoTrans=111, CblasTrans=112,
                             CblasConjTrans=113, CblasConjNoTrans=114};

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CBLAS_ALPHA;
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef const CBLAS_FLOAT *         CBLAS_ALPHA;
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               CBLAS_INT m, CBLAS_INT n,
               CBLAS_ALPHA _alpha,
               const CBLAS_FLOAT  *_x, CBLAS_INT incX,
               const CBLAS_FLOAT  *_y, CBLAS_INT incY,
               CBLAS_FLOAT  *_A, CBLAS_INT ldA)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;


        const CXXBLAS_FLOAT *x = reinterpret_cast<const CXXBLAS_FLOAT *>(_x);
        const CXXBLAS_FLOAT *y = reinterpret_cast<const CXXBLAS_FLOAT *>(_y);
        CXXBLAS_FLOAT *A = reinterpret_cast<CXXBLAS_FLOAT *>(_A);

    #ifndef COMPLEX_FLOAT1
        CXXBLAS_FLOAT alpha = _alpha;
    #else
        CXXBLAS_FLOAT alpha(_alpha[0], _alpha[1]);
    #endif

        CXXBLAS_FLOAT *__A = A;
        CBLAS_INT      __ldA = ldA;
        StorageOrder   __order = order;

#   ifdef TEST_ROW_MAJOR
        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocateFullStorage(__order, m, n, __A, __ldA);
        switchFullStorageOrder(order, m, n, A, ldA, __A, __ldA);
#   endif

        cxxblas::geru(__order,
                      m, n,
                      alpha,
                      x, incX,
                      y, incY,
                      __A, __ldA);

#   ifdef TEST_ROW_MAJOR
        switchFullStorageOrder(__order, m, n, __A, __ldA, A, ldA);
        releaseStorage(__A);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const CBLAS_INT *_m, const CBLAS_INT *_n,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *x, const CBLAS_INT *_incX,
              const CBLAS_FLOAT *y, const CBLAS_INT *_incY,
              CBLAS_FLOAT *A, const CBLAS_INT *_ldA)
    {
        CBLAS_INT m       = *_m;
        CBLAS_INT n       = *_n;
        CBLAS_INT incX    = *_incX;
        CBLAS_INT incY    = *_incY;
        CBLAS_INT ldA     = *_ldA;

#   ifndef COMPLEX_FLOAT1
        CBLAS_ALPHA alpha  = *_alpha;
#   else
        CBLAS_ALPHA alpha = _alpha;
#   endif

        CBLAS_INT info = 0;
        if (incY==0) {
            info = 7;
        }
        if (incX==0) {
            info = 5; 
        }
        if (ldA<std::max(1, m)) {
            info = 9;
        }
        if (n<0) {
            info = 2;
        }
        if (m<0) {
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
        CBLAS_NAME(CblasColMajor,
                   m, n,
                   alpha,
                   x, incX,
                   y, incY,
                   A, ldA);
    }
#endif // CREATE_BLAS

} // extern "C"