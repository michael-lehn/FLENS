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

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               CBLAS_INT m, CBLAS_INT n,
               CBLAS_FLOAT alpha,
               const CBLAS_FLOAT  *x, CBLAS_INT incX,
               const CBLAS_FLOAT  *y, CBLAS_INT incY,
               CBLAS_FLOAT  *A, CBLAS_INT ldA)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;

        CBLAS_FLOAT *__A = A;
        CBLAS_INT    __ldA = ldA;
        StorageOrder __order = order;

#   ifdef TEST_ROW_MAJOR
        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocateFullStorage(__order, m, n, __A, __ldA);
        switchFullStorageOrder(order, m, n, A, ldA, __A, __ldA);
#   endif

        cxxblas::ger(__order,
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
        CBLAS_FLOAT alpha = *_alpha;
        CBLAS_INT incX    = *_incX;
        CBLAS_INT incY    = *_incY;
        CBLAS_INT ldA     = *_ldA;

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