#include <cxxblas/cxxblas.cxx>
#include <cxxblas/interface/aux.h>

using cxxblas::StorageOrder;
using cxxblas::ColMajor;
using cxxblas::RowMajor;

using cxxblas::StorageUpLo;
using cxxblas::Upper;
using cxxblas::Lower;

extern "C" {

    void xerbla_(const char* srname, int* info);

#ifndef COMPLEX_FLOAT1
    typedef CBLAS_FLOAT                 CXXBLAS_FLOAT;
#else
    typedef std::complex<CBLAS_FLOAT>   CXXBLAS_FLOAT;
#endif

    enum CBLAS_ORDER        {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_UPLO         {CblasUpper=121, CblasLower=122};

    void
    CBLAS_NAME(enum CBLAS_ORDER _order,
               enum CBLAS_UPLO _upLo,
               CBLAS_INT n,
               CBLAS_FLOAT alpha,
               const CBLAS_FLOAT  *_x, CBLAS_INT incX,
               CBLAS_FLOAT  *_A, CBLAS_INT ldA)
#ifdef CREATE_CBLAS
    {
        StorageOrder order = (_order==CblasColMajor) ? ColMajor
                                                     : RowMajor;
        StorageUpLo upLo = Lower;
        if (_upLo==CblasUpper) {
            upLo = Upper;
        }

        const CXXBLAS_FLOAT *x = reinterpret_cast<const CXXBLAS_FLOAT *>(_x);
        CXXBLAS_FLOAT *A = reinterpret_cast<CXXBLAS_FLOAT *>(_A);

        CXXBLAS_FLOAT *__A = A;
        CBLAS_INT      __ldA = ldA;
        StorageOrder   __order = order;

#   ifdef TEST_ROW_MAJOR
        __order = (order==ColMajor) ? RowMajor : ColMajor;
        allocateFullStorage(__order, n, n, __A, __ldA);
        switchFullStorageOrder(order, n, n, A, ldA, __A, __ldA);
#   endif

        cxxblas::her(__order, upLo,
                     n,
                     alpha,
                     x, incX,
                     __A, __ldA);

#   ifdef TEST_ROW_MAJOR
        switchFullStorageOrder(__order, n, n, __A, __ldA, A, ldA);
        releaseStorage(__A);
#   endif
    }
#else
    ;
#endif // CREATE_CBLAS

#ifdef CREATE_BLAS
    void
    BLAS_NAME(const char *_upLo,
              const CBLAS_INT *_n,
              const CBLAS_FLOAT *_alpha,
              const CBLAS_FLOAT *x, const CBLAS_INT *_incX,
              CBLAS_FLOAT *A, const CBLAS_INT *_ldA)
    {
        bool checkUpLo = false;
        CBLAS_UPLO upLo = CblasUpper;
        if ((*_upLo=='L') || (*_upLo=='l')) {
            upLo = CblasLower;
            checkUpLo = true;
        }
        if ((*_upLo=='U') || (*_upLo=='u')) {
            checkUpLo = true;
        }

        CBLAS_INT n       = *_n;
        CBLAS_FLOAT alpha = *_alpha;
        CBLAS_INT incX    = *_incX;
        CBLAS_INT ldA     = *_ldA;


        CBLAS_INT info = 0;
        if (incX==0) {
            info = 5; 
        }
        if (ldA<std::max(1, n)) {
            info = 7;
        }
        if (n<0) {
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
        CBLAS_NAME(CblasColMajor, upLo,
                   n,
                   alpha,
                   x, incX,
                   A, ldA);
    }
#endif // CREATE_BLAS

} // extern "C"